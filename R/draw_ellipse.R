#' Taken from plotrix
.draw_ellipse <-
function (x, y, a = 1, b = 1, angle = 0, segment = NULL, arc.only = TRUE,
    deg = TRUE, nv = 100, border = NULL, col = NA, lty = 1, lwd = 1,
    ...)
{
    if (is.null(segment)) {
        if (deg) {
            segment <- c(0, 360)
        } else {
            segment <- c(0, 2 * pi)
        }
    }
    draw1ellipse <- function(x, y, a = 1, b = 1, angle = 0, segment = NULL,
        arc.only = TRUE, nv = 100, deg = TRUE, border = NULL,
        col = NA, lty = 1, lwd = 1, ...) {
        if (deg) {
            angle <- angle * pi/180
            segment <- segment * pi/180
        }
        z <- seq(segment[1], segment[2], length = nv + 1)
        xx <- a * cos(z)
        yy <- b * sin(z)
        alpha <- xyangle(xx, yy, directed = TRUE, deg = FALSE)
        rad <- sqrt(xx^2 + yy^2)
        xp <- rad * cos(alpha + angle) + x
        yp <- rad * sin(alpha + angle) + y
        if (!arc.only) {
            xp <- c(x, xp, x)
            yp <- c(y, yp, y)
        }
        polygon(xp, yp, border = border, col = col, lty = lty,
            lwd = lwd, ...)
        invisible(NULL)
    }
    xyangle <- function(x, y, directed = FALSE, deg = TRUE) {
        if (missing(y)) {
            y <- x[, 2]
            x <- x[, 1]
        }
        out <- atan2(y, x)
        if (!directed)
            out <- out%%pi
        if (deg)
            out <- out * 180/pi
        out
    }
    if (missing(y)) {
        y <- x[, 2]
        x <- x[, 1]
    }
    n <- length(x)
    if (length(a) < n)
        a <- rep(a, n)[1:n]
    if (length(b) < n)
        b <- rep(b, n)[1:n]
    if (length(angle) < n)
        angle <- rep(angle, n)[1:n]
    if (length(col) < n)
        col <- rep(col, n)[1:n]
    if (length(border) < n)
        border <- rep(border, n)[1:n]
    if (length(nv) < n)
        nv <- rep(nv, n)[1:n]
    if (n == 1) {
        draw1ellipse(x, y, a, b, angle = angle, segment = segment,
            arc.only = arc.only, deg = deg, nv = nv, col = col,
            border = border, lty = lty, lwd = lwd, ...)
    } else {
        if (length(segment) < 2 * n)
            segment <- matrix(rep(segment, n), n, 2, byrow = TRUE)
        lapply(1:n, function(i) draw1ellipse(x[i], y[i], a[i],
            b[i], angle = angle[i], segment = segment[i, ], arc.only = arc.only,
            deg = deg, nv = nv[i], col = col[i], border = border[i],
            lty = lty, lwd = lwd, ...))
    }
    invisible(NULL)
}
