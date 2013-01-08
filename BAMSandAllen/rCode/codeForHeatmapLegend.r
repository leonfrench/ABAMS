    ##------ Some palettes ------------
demo.pal <-
       function(n, border = if (n<32) "light gray" else NA,
                main = paste("color palettes;  n=",n),
                ch.col = c( "heat.colors(n)"))
     {
         nt <- length(ch.col)
         i <- 1:n; j <- n / nt; d <- j/6; dy <- 2*d
         plot(i,i+d, type="n", yaxt="n", ylab="", main=main)
         for (k in 1:nt) {
             rect(i-.5, (k-1)*j+ dy, i+.5, k*j,
                  col = eval(parse(text=ch.col[k])), border = border)
             text(2*j,  k * j +dy/4, ch.col[k])
         }
     }
     demo.pal(64)

