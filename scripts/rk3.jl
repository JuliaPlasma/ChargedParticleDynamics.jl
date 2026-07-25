
function RK3()
    a = [[1 / 2 -√3 / 6]
        [+√3 / 6 1 / 2]]
    b = [1 / 2, 1 / 2]
    c = [1 / 2, 1 / 2]
    o = 3

    RK(Tableau{Float64}(:RK3, o, a, b, c))
end
