#### Fit light curves to multiperiodic pulsators 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(glmnet)
library(Hmisc)

phases <- seq(0, 0.99, 0.01)
nlambda <- 1000
N_max <- 8

Fourier <- function(t=phases, p=1, Nmax=N_max) {
    base_freqs <- 1/p
    combos <- expand.grid(rep(list(-Nmax:Nmax), length(base_freqs)))
    freqs <- apply(combos, 1, function(combo) { base_freqs %*% combo })
    freqs <- unique(abs(freqs))
    freqs <- rev(freqs[freqs>0])
    num_freqs <- length(freqs)
    X <- matrix(nrow=length(t), ncol=2*num_freqs) 
    for (ii in 1:num_freqs) {
        X[,2*(ii-1)+1] <- sin(2*pi*t*freqs[ii])
        X[,2*(ii-1)+2] <- cos(2*pi*t*freqs[ii])
    }
    return(list(X=X, freqs=freqs))
}
evenly_spaced <- Fourier()$X

start_devices <- function(save_dir, filename, plot_name, output_fmt) {
    dir.create(save_dir, showWarnings=FALSE)
    id <- sub('.+/', '', sub('\\..+', '', filename))
    if (grepl('pdf', output_fmt)) {
        cairo_pdf(file.path(save_dir, paste0(id, '_', plot_name, '.pdf')), 
                  width=6, height=6, family="Palatino")
    } else if (grepl('png', output_fmt)) {
        png(file.path(save_dir, paste0(id, '_', plot_name, '.png')),
            width=600, height=800, family="Palatino")
        par(cex=1.25)
    }
}

plot_residuals <- function(x, m_res, err, xlab="Phase") {
    errbar(x, m_res, 
           m_res+err/2, 
           m_res-err/2, 
           xlim=round(range(x)), xaxs='i', tck=0.01, cex=0, cap=0,
           main='', ylab=expression(m - hat(m)), xlab=xlab)
    title(main='Residual plot')
    lines(c(round(max(x)), 0), c(0, 0), col='red', lty=2)
    minor.tick(nx=5, ny=5, tick.ratio=-0.25)
}

plot_single <- function(period, photometry, phase, m_even, m_hat, filename) {
    ## Phase plot
    par(las=1, mar=c(4, 4, 3, 1), mgp=c(2.75, 0.25, 0), mfrow=c(3, 1), 
        cex=par()$cex)
    errbar(phase, photometry$m, 
           photometry$m+photometry$e/2, 
           photometry$m-photometry$e/2,
           xlim=c(0,1), ylim=rev(range(photometry$m, m_even)),
           xaxs='i', tck=0.01, cex=0, cap=0,
           main='', ylab='Magnitude', xlab="Phase")
    title(main=paste0(sub('\\..+', '', filename), 
                      " (", period, " day period)"))
    lines(c(phases, 1), c(m_even, m_even[1]), col='red')
    minor.tick(nx=5, ny=5, tick.ratio=-0.25)
    
    ## Residual plot
    m_res <- photometry$m-m_hat
    plot_residuals(phase, m_res, photometry$e)
    
    ## Normal QQ plot
    qqnorm(m_res, pch=3)
    qqline(m_res)
}

plot_multiple <- function(periods, photometry, m_hat, filename, t0s=0) {
    ## phase plots
    layout(matrix(c(if (length(periods) %% 2 == 0) 1 else c(1,1),
                    2:length(periods), length(periods)+1, length(periods)+2),
                  ncol=2, byrow=TRUE))
    par(las=1, mar=c(4, 4, 3, 1), mgp=c(2.75, 0.25, 0), cex=par()$cex)
    for (ii in 1:length(periods)) {
        t0 <- ifelse(length(t0s) == length(periods), t0s[ii], t0)
        phase <- ((photometry$t-t0) %% periods[ii])/periods[ii]
        errbar(phase, photometry$m, 
               photometry$m+photometry$e/2, 
               photometry$m-photometry$e/2,
               xlim=c(0,1), ylim=rev(range(photometry$m)),
               xaxs='i', tck=0.01, cex=0, cap=0,
               main='', ylab='Magnitude', 
               xlab=paste0("Phase (", periods[ii], " day period)"))
        if (ii==1) title(main=sub('\\..+', '', filename))
        minor.tick(nx=5, ny=5, tick.ratio=-0.25)
    }
    
    ## Residual plot
    m_res <- photometry$m-m_hat
    plot_residuals(1:length(m_res), m_res, photometry$e, "Index")
    
    ## Normal QQ plot
    qqnorm(m_res, pch=3)
    qqline(m_res)
}

plot_time_series <- function(photometry, periods, m_hat, cvfit) {
    # find regions where there aren't gaps 
    # plotting this is far more of an art than a science... 
    gaps <- diff(photometry$t)
    big_gaps <- unique(c(0, which(gaps > median(gaps) + 10*mad(gaps)), 
        length(photometry$t)))
    should_plot <- diff(big_gaps) > median(diff(big_gaps))
    num_plots <- sum(should_plot)
    par(las=1, mar=c(4, 4, 3, 1), mgp=c(2.75, 0.25, 0), 
        mfrow=c(round((num_plots/2)+0.49), 2), cex=par()$cex)
    for (ii in 1:(length(big_gaps)-1)) {
        if (!should_plot[ii]) next
        start_i <- big_gaps[ii]+1
        end_i <- big_gaps[ii+1]
        errbar(photometry$t[start_i:end_i], photometry$m[start_i:end_i], 
               photometry$m[start_i:end_i]+photometry$e[start_i:end_i]/2, 
               photometry$m[start_i:end_i]-photometry$e[start_i:end_i]/2,
               ylim=rev(range(photometry$m, m_hat)),
               tck=0.01, cex=0, cap=0,
               main='', ylab='Magnitude', 
               xlab="HJD-2450000")
        minor.tick(nx=5, ny=5, tick.ratio=-0.25)
        t_even <- seq(photometry$t[start_i], photometry$t[end_i], 0.01)
        Fourier_even <- Fourier(t_even, periods)
        m_even <- predict(cvfit, newx=Fourier_even$X, 
            s="lambda.min", exact=TRUE)
        lines(t_even, m_even, col=rgb(1,0,0,0.5))
    }
}

fit_single <- function(photometry, period=1, t0=0, n_lambda=nlambda) {
    phase <- ((photometry$t-t0) %% period)/period
    Fourier_phase <- Fourier(phase, 1)$X
    cvfit <- cv.glmnet(Fourier_phase, photometry$m, nlambda=n_lambda)
    m_hat <- predict(cvfit, newx=Fourier_phase, s="lambda.min", exact=TRUE)
    m_even <- predict(cvfit, newx=evenly_spaced, s="lambda.min", exact=TRUE)
    return(list(m_even=m_even, m_hat=m_hat, cvfit=cvfit, phase=phase))
}

fit_multiple <- function(photometry, periods, n_lambda=nlambda) {
    Fourier_space <- Fourier(photometry$t, periods)
    cvfit <- cv.glmnet(Fourier_space$X, photometry$m, nlambda=n_lambda)
    m_hat <- predict(cvfit, newx=Fourier_space$X, s="lambda.min", exact=TRUE)
    return(list(m_hat=m_hat, cvfit=cvfit, freqs=Fourier_space$freqs,
                coefs=coef(cvfit)))
}

fit_lightcurve <- function(filename, periods=1, t0s=0, show_plot=TRUE, 
        save_dir=NULL, output_fmt='.png', n_lambda=nlambda, 
        plot_ts=TRUE) {
    if (!file.exists(filename)) 
        return("Cannot find file")
    photometry <- read.table(filename, col.names=c('t', 'm', 'e'))
    
    if (show_plot && !is.null(save_dir))
        start_devices(save_dir, filename, 'phase', output_fmt)
    
    if (length(periods) == 1) {
        fit <- fit_single(photometry, periods, t0=t0s[1], n_lambda)
        if (show_plot) 
            plot_single(periods, photometry, fit$phase, 
                        fit$m_even, fit$m_hat, filename)
    } else { 
        fit <- fit_multiple(photometry, periods, n_lambda)
        if (show_plot) 
            plot_multiple(periods, photometry, fit$m_hat, filename, t0s)
    }
    
    if (show_plot) {
        if (!is.null(save_dir)) dev.off()
        if (plot_ts) {
            if (!is.null(save_dir))
                start_devices(save_dir, filename, 'ts', output_fmt)
            else dev.new()
            plot_time_series(photometry, periods, fit$m_hat, fit$cvfit)
            if (!is.null(save_dir)) dev.off()
        }
    }
    return(fit)
}

fit_lightcurve('OGLE-LMC-CEP-0002.dat', 3.118120, 2171.239,
               save_dir='multiplots', output_fmt='.png')

fit_lightcurve('OGLE-SMC-CEP-0408.dat', 
               c(1.7901765, 1.3140538), c(624.75372, 624.97043),
               save_dir='multiplots', output_fmt='.png')

fit_lightcurve('OGLE-SMC-CEP-3867.dat', 
               c(0.2688471, 0.2173800, 0.1824204), 
               c(2104.60491, 2104.81172, 2104.72340),
               save_dir='multiplots', output_fmt='.png')
