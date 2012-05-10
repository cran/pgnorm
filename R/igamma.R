igamma <-
function(x, a)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Incomplete Gamma Function P(a, x) with 
    #   Re(a) > 0 for complex or real argument "x" and for 
    #   complex or real index "z" 
    
    # Arguments:
    #   z - a complex or real vector
    #   a - a complex or real numeric value
    
    # Details:
    #   [AS] formula 6.5.1  
    #   $ frac{1}{Gamma(a)}  * \int_0^x e^{-t} t^{a-1} dt $
    
    # FUNCTION:
    
    # igamma:
    if (!is.complex(x) && !is.complex(a)) {
        # Use R's pgamma() function:
        # if (a < 0) Not suppported ...
        result = pgamma(x, a)  
    } else {
        # Why not derive the result from KummersM ?
        log = FALSE
        if (log) {
            # Not yet supported:
            result = kummerM(a, a + 1, -x, lnchf = 1) + a*log(x) - log(a) 
        } else {
            result = kummerM(a, a + 1, -x, lnchf = 0) * x^a / a 
        } 
    }
        
    # Return Value:
    result
}
