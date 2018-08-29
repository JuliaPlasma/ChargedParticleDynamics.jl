
using Printf

function power10ticks(x, xmax, pos)
    if x == 0
        return "\$ 0 \$"
    end

    xpower      = log10(xmax)
    exponent    = @sprintf("%2d",   floor(Int64, xpower))
    coefficient = @sprintf("%2.1f", x / (10^floor(xpower)))
    ticks       = "\$ $coefficient \\times 10^{ $exponent } \$"

    return ticks
end


function get_default_formatter(t)

    xt(x, pos) = power10ticks(x, t[end], pos)

    xf = matplotlib[:ticker][:FuncFormatter](xt)

    yf = matplotlib[:ticker][:ScalarFormatter]()
    yf[:set_powerlimits]((-1,+1))
    yf[:set_scientific](true)
    yf[:set_useOffset](true)

    return (xf, yf)
end
