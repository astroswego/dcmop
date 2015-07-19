#### Fit light curves to multiperiodic pulsators 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(glmnet)
library(Hmisc)

phases <- seq(0, 0.99, 0.01)
n_lambda <- 1000
N_max <- 8

Fourier <- function(t=phases, p=1, Nmax=N_max) {
    f <- 1/p
    combos <- expand.grid(rep(list(-Nmax:Nmax), length(f)))
    combo_freqs <- apply(combos, 1, function(combo) { f %*% combo })
    combo_freqs <- unique(abs(combo_freqs))
    combo_freqs <- rev(combo_freqs[combo_freqs>0])
    num_freqs <- length(combo_freqs)
    X <- matrix(nrow=length(t), ncol=2*num_freqs) 
    for (ii in 1:num_freqs) {
        X[,2*(ii-1)+1] <- sin(2*pi*t*combo_freqs[ii])
        X[,2*(ii-1)+2] <- cos(2*pi*t*combo_freqs[ii])
    }
    return(list(X=X, combo_freqs=combo_freqs))
}
evenly_spaced <- Fourier()$X

fit_single_lightcurve <- function(filename, period=1, t0=0, show_plot=TRUE,
        n_lambda.=n_lambda) {
    if (!file.exists(filename)) return(rep(NA, length(phases)))
    photometry <- read.table(filename, col.names=c('t', 'm', 'e'))
    phase <- ((photometry$t-t0) %% period)/period
    Fourier_phase <- Fourier(phase, 1)$X
    cvfit <- cv.glmnet(Fourier_phase, photometry$m, nlambda=n_lambda.)
    m_hat <- predict(cvfit, newx=Fourier_phase, s="lambda.min", exact=TRUE)
    m_even <- predict(cvfit, newx=evenly_spaced, s="lambda.min", exact=TRUE)
    if (show_plot) {
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
        errbar(phase, m_res, 
               m_res+photometry$e/2, 
               m_res-photometry$e/2, 
               xlim=c(0,1), xaxs='i', tck=0.01, cex=0, cap=0,
               main='', ylab=expression(m - hat(m)), xlab="Phase")
        title(main='Residual plot')
        lines(c(1, 0), c(0, 0), col='black', lty=2)
        minor.tick(nx=5, ny=5, tick.ratio=-0.25)
        
        ## Normal QQ plot
        qqnorm(m_res, pch=3)
        qqline(m_res)
    }
    return(m_even)
}

fit_lightcurve <- function(filename, periods=1, t0s=0, 
        show_plot=TRUE, save_dir=NULL, output='png', n_lambda.=n_lambda) {
    if (!file.exists(filename)) return("Cannot find file")
    
    ## Prepare plotting devices 
    if (show_plot && !is.null(save_dir)) {
        dir.create(save_dir, showWarnings=FALSE)
        id <- sub('.+/', '', sub('\\..+', '', filename))
        if (grepl('pdf', output)) {
            cairo_pdf(file.path(save_dir, paste0(id, '.pdf')), 
                      width=6, height=6, family="Palatino")
        } else if (grepl('png', output)) {
            png(file.path(save_dir, paste0(id, '.png')),
                width=600, height=800, family="Palatino")
            par(cex=1.25)
        }
    }
    
    ## if there's only one period, fit it and exit
    if (length(periods) == 1) {
        m_even <- fit_single_lightcurve(filename, periods, t0s[1], show_plot,
            n_lambda.)
        if (show_plot && !is.null(save_dir)) dev.off()
        return(m_even)
    } 
    
    # otherwise, fit time series and make phase plots for each period 
    photometry <- read.table(filename, col.names=c('t', 'm', 'e'))
    Fourier_space <- Fourier(photometry$t, periods)
    cvfit <- cv.glmnet(Fourier_space$X, photometry$m, nlambda=n_lambda.)
    m_hat <- predict(cvfit, newx=Fourier_space$X, s="lambda.min", exact=TRUE)
    
    if (show_plot) {
        ## phase plots
        if (length(periods) %% 2 == 0) {
            #par(mfcol=c((length(periods)+2)/2, 2), cex=par()$cex)
            layout(matrix(c(1:length(periods), 
                            length(periods)+1, length(periods)+2),
                          ncol=2, byrow=TRUE))
        } else {
            layout(matrix(c(rep(1:length(periods), each=2), 
                            rep(length(periods)+1, 3),
                            rep(length(periods)+2, 3)),
                          ncol=6, byrow=TRUE))
        }
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
        errbar(1:length(m_res), m_res, 
               m_res+photometry$e/2, 
               m_res-photometry$e/2, 
               xaxs='i', tck=0.01, cex=0, cap=0,
               main='', ylab=expression(m - hat(m)), xlab="Index")
        title(main='Residual plot')
        lines(c(length(m_res), 0), c(0, 0), col='black', lty=2)
        minor.tick(nx=5, ny=5, tick.ratio=-0.25)
        
        ## Normal QQ plot
        qqnorm(m_res, pch=3)
        qqline(m_res)
        
        ## prepare devices again to plot time series 
        if (!is.null(save_dir)) {
            dev.off()
            if (grepl('pdf', output)) {
                cairo_pdf(file.path(save_dir, paste0(id, '_ts.pdf')), 
                          width=12, height=12, family="Palatino")
            } else if (grepl('png', output)) {
                png(file.path(save_dir, paste0(id, '_ts.png')),
                    width=1200, height=1600, family="Palatino")
                par(cex=1.25)
            }
        } else {
            dev.new()
        }
        
        ## plot time series and residuals 
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
                   xaxs='i', tck=0.01, cex=0, cap=0,
                   main='', ylab='Magnitude', 
                   xlab="HJD-2450000")
            minor.tick(nx=5, ny=5, tick.ratio=-0.25)
            t_even <- seq(photometry$t[start_i], photometry$t[end_i], 0.01)
            Fourier_even <- Fourier(t_even, periods)
            m_even <- predict(cvfit, newx=Fourier_even$X, 
                s="lambda.min", exact=TRUE)
            lines(t_even, m_even, col=rgb(1,0,0,0.5))
        }
        if (!is.null(save_dir)) dev.off()
    }
}

#fit_lightcurve('OGLE-LMC-CEP-0002.dat', 3.118120, 2171.239,
#               save_dir='multiplots', output='.png')

fit_lightcurve('OGLE-SMC-CEP-0408.dat', 
               c(1.7901765, 1.3140538), c(624.75372, 624.97043),
               save_dir='multiplots', output='.png')

fit_lightcurve('OGLE-SMC-CEP-3867.dat', 
               c(0.2688471, 0.2173800, 0.1824204), 
               c(2104.60491, 2104.81172, 2104.72340),
               save_dir='multiplots', output='.pdf')

