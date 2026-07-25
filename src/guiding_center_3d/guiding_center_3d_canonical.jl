
d²v₁dq₁dq₁(t, q, p) = -d²A₁dx₁dx₁(t, q)
d²v₁dq₁dq₂(t, q, p) = -d²A₁dx₁dx₂(t, q)
d²v₁dq₁dq₃(t, q, p) = -d²A₁dx₁dx₃(t, q)
d²v₁dq₂dq₂(t, q, p) = -d²A₁dx₂dx₂(t, q)
d²v₁dq₂dq₃(t, q, p) = -d²A₁dx₂dx₃(t, q)
d²v₁dq₃dq₃(t, q, p) = -d²A₁dx₃dx₃(t, q)

d²v₂dq₁dq₁(t, q, p) = -d²A₂dx₁dx₁(t, q)
d²v₂dq₁dq₂(t, q, p) = -d²A₂dx₁dx₂(t, q)
d²v₂dq₁dq₃(t, q, p) = -d²A₂dx₁dx₃(t, q)
d²v₂dq₂dq₂(t, q, p) = -d²A₂dx₂dx₂(t, q)
d²v₂dq₂dq₃(t, q, p) = -d²A₂dx₂dx₃(t, q)
d²v₂dq₃dq₃(t, q, p) = -d²A₂dx₃dx₃(t, q)

d²v₃dq₁dq₁(t, q, p) = -d²A₃dx₁dx₁(t, q)
d²v₃dq₁dq₂(t, q, p) = -d²A₃dx₁dx₂(t, q)
d²v₃dq₁dq₃(t, q, p) = -d²A₃dx₁dx₃(t, q)
d²v₃dq₂dq₂(t, q, p) = -d²A₃dx₂dx₂(t, q)
d²v₃dq₂dq₃(t, q, p) = -d²A₃dx₂dx₃(t, q)
d²v₃dq₃dq₃(t, q, p) = -d²A₃dx₃dx₃(t, q)


d²Hdq₁dq₁(t, q, p, params) = (
    dv₁dq₁(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dq₁(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dq₁(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₁dq₁(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₁dq₁(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₁dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dq₁(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dq₁(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dq₁(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₁dx₁(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₁dx₁(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₁dx₁(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₁dx₁(t, q) -
    dE₁dx₁(t, q)
)

d²Hdq₁dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₁dq₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₁dq₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₁dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dq₂(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dq₂(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dq₂(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₁dx₂(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₁dx₂(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₁dx₂(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₁dx₂(t, q) -
    dE₁dx₂(t, q)
)

d²Hdq₁dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₁dq₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₁dq₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₁dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dq₃(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dq₃(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dq₃(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₁dx₃(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₁dx₃(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₁dx₃(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₁dx₃(t, q) -
    dE₁dx₃(t, q)
)

d²Hdq₁dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₁(t, q, p)
)

d²Hdq₁dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₂(t, q, p)
)

d²Hdq₁dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dq₁(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dq₁(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₃(t, q, p)
)


d²Hdq₂dq₁(t, q, p, params) = (
    dv₁dq₁(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dq₁(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dq₁(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₁dq₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₁dq₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₁dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dq₁(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dq₁(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dq₁(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₂dx₁(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₂dx₁(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₂dx₁(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₂dx₁(t, q) -
    dE₂dx₁(t, q)
)

d²Hdq₂dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₂dq₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₂dq₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₂dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dq₂(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dq₂(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dq₂(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₂dx₂(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₂dx₂(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₂dx₂(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₂dx₂(t, q) -
    dE₂dx₂(t, q)
)

d²Hdq₂dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₂dq₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₂dq₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₂dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dq₃(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dq₃(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dq₃(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₂dx₃(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₂dx₃(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₂dx₃(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₂dx₃(t, q) -
    dE₂dx₃(t, q)
)

d²Hdq₂dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₁(t, q, p)
)

d²Hdq₂dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₂(t, q, p)
)

d²Hdq₂dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dq₂(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dq₂(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₃(t, q, p)
)


d²Hdq₃dq₁(t, q, p, params) = (
    dv₁dq₁(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dq₁(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dq₁(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₃dq₁(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₃dq₁(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₃dq₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dq₁(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dq₁(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dq₁(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₃dx₁(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₃dx₁(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₃dx₁(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₃dx₁(t, q) -
    dE₃dx₁(t, q)
)

d²Hdq₃dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₃dq₂(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₃dq₂(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₃dq₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dq₂(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dq₂(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dq₂(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₃dx₂(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₃dx₂(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₃dx₂(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₃dx₂(t, q) -
    dE₃dx₂(t, q)
)

d²Hdq₃dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * g¹¹(t, q) * d²v₁dq₃dq₃(t, q, p) +
    v₂(t, q, p) * g²²(t, q) * d²v₂dq₃dq₃(t, q, p) +
    v₃(t, q, p) * g³³(t, q) * d²v₃dq₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dq₃(t, q, p) * 2 +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dq₃(t, q, p) * 2 +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dq₃(t, q, p) * 2 +
    v₁(t, q, p) * d²g¹¹dx₃dx₃(t, q) * v₁(t, q, p) / 2 +
    v₂(t, q, p) * d²g²²dx₃dx₃(t, q) * v₂(t, q, p) / 2 +
    v₃(t, q, p) * d²g³³dx₃dx₃(t, q) * v₃(t, q, p) / 2 +
    params.μ * d²Bdx₃dx₃(t, q) -
    dE₃dx₃(t, q)
)

d²Hdq₃dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₁(t, q, p)
)

d²Hdq₃dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₂(t, q, p)
)

d²Hdq₃dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dq₃(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dq₃(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dq₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₃(t, q, p)
)


d²Hdp₁dq₁(t, q, p, params) = (
    dv₁dq₁(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dq₁(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dq₁(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dq₁(t, q, p, params) = (
    dv₁dq₁(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dq₁(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dq₁(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₂(t, q, p)
)

d²Hdp₃dq₁(t, q, p, params) = (
    dv₁dq₁(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    dv₂dq₁(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    dv₃dq₁(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₁(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₁(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₁(t, q) * dv₃dp₃(t, q, p)
)

d²Hdp₁dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)
)

d²Hdp₃dp₁(t, q, p, params) = (
    dv₁dp₁(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    dv₂dp₁(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    dv₃dp₁(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p)
)


d²Hdp₁dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₂(t, q, p)
)

d²Hdp₃dq₂(t, q, p, params) = (
    dv₁dq₂(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    dv₂dq₂(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    dv₃dq₂(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₂(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₂(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₂(t, q) * dv₃dp₃(t, q, p)
)

d²Hdp₁dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)
)

d²Hdp₃dp₂(t, q, p, params) = (
    dv₁dp₂(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    dv₂dp₂(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    dv₃dp₂(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p)
)


d²Hdp₁dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₁(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₁(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₂(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₂(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₂(t, q, p)
)

d²Hdp₃dq₃(t, q, p, params) = (
    dv₁dq₃(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    dv₂dq₃(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    dv₃dq₃(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p) +
    v₁(t, q, p) * dg¹¹dx₃(t, q) * dv₁dp₃(t, q, p) +
    v₂(t, q, p) * dg²²dx₃(t, q) * dv₂dp₃(t, q, p) +
    v₃(t, q, p) * dg³³dx₃(t, q) * dv₃dp₃(t, q, p)
)

d²Hdp₁dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dp₁(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dp₁(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dp₁(t, q, p)
)

d²Hdp₂dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dp₂(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dp₂(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dp₂(t, q, p)
)

d²Hdp₃dp₃(t, q, p, params) = (
    dv₁dp₃(t, q, p) * g¹¹(t, q) * dv₁dp₃(t, q, p) +
    dv₂dp₃(t, q, p) * g²²(t, q) * dv₂dp₃(t, q, p) +
    dv₃dp₃(t, q, p) * g³³(t, q) * dv₃dp₃(t, q, p)
)



d²g₁dq₁dq₁(t, q, p) = (d²b₁dx₁dx₁(t, q) * v₂(t, q, p) + db₁dx₁(t, q) * dv₂dq₁(t, q, p)) -
                      (d²b₂dx₁dx₁(t, q) * v₁(t, q, p) + db₂dx₁(t, q) * dv₁dq₁(t, q, p)) +
                      (db₁dx₁(t, q) * dv₂dq₁(t, q, p) + b₁(t, q) * d²v₂dq₁dq₁(t, q, p)) -
                      (db₂dx₁(t, q) * dv₁dq₁(t, q, p) + b₂(t, q) * d²v₁dq₁dq₁(t, q, p))

d²g₁dq₁dq₂(t, q, p) = (d²b₁dx₁dx₂(t, q) * v₂(t, q, p) + db₁dx₂(t, q) * dv₂dq₁(t, q, p)) -
                      (d²b₂dx₁dx₂(t, q) * v₁(t, q, p) + db₂dx₂(t, q) * dv₁dq₁(t, q, p)) +
                      (db₁dx₁(t, q) * dv₂dq₂(t, q, p) + b₁(t, q) * d²v₂dq₁dq₂(t, q, p)) -
                      (db₂dx₁(t, q) * dv₁dq₂(t, q, p) + b₂(t, q) * d²v₁dq₁dq₂(t, q, p))

d²g₁dq₁dq₃(t, q, p) = (d²b₁dx₁dx₃(t, q) * v₂(t, q, p) + db₁dx₃(t, q) * dv₂dq₁(t, q, p)) -
                      (d²b₂dx₁dx₃(t, q) * v₁(t, q, p) + db₂dx₃(t, q) * dv₁dq₁(t, q, p)) +
                      (db₁dx₁(t, q) * dv₂dq₃(t, q, p) + b₁(t, q) * d²v₂dq₁dq₃(t, q, p)) -
                      (db₂dx₁(t, q) * dv₁dq₃(t, q, p) + b₂(t, q) * d²v₁dq₁dq₃(t, q, p))

d²g₁dq₂dq₁(t, q, p) = (d²b₁dx₂dx₁(t, q) * v₂(t, q, p) + db₁dx₁(t, q) * dv₂dq₂(t, q, p)) -
                      (d²b₂dx₂dx₁(t, q) * v₁(t, q, p) + db₂dx₁(t, q) * dv₁dq₂(t, q, p)) +
                      (db₁dx₂(t, q) * dv₂dq₁(t, q, p) + b₁(t, q) * d²v₂dq₁dq₂(t, q, p)) -
                      (db₂dx₂(t, q) * dv₁dq₁(t, q, p) + b₂(t, q) * d²v₁dq₁dq₂(t, q, p))

d²g₁dq₂dq₂(t, q, p) = (d²b₁dx₂dx₂(t, q) * v₂(t, q, p) + db₁dx₂(t, q) * dv₂dq₂(t, q, p)) -
                      (d²b₂dx₂dx₂(t, q) * v₁(t, q, p) + db₂dx₂(t, q) * dv₁dq₂(t, q, p)) +
                      (db₁dx₂(t, q) * dv₂dq₂(t, q, p) + b₁(t, q) * d²v₂dq₂dq₂(t, q, p)) -
                      (db₂dx₂(t, q) * dv₁dq₂(t, q, p) + b₂(t, q) * d²v₁dq₂dq₂(t, q, p))

d²g₁dq₂dq₃(t, q, p) = (d²b₁dx₂dx₃(t, q) * v₂(t, q, p) + db₁dx₃(t, q) * dv₂dq₂(t, q, p)) -
                      (d²b₂dx₂dx₃(t, q) * v₁(t, q, p) + db₂dx₃(t, q) * dv₁dq₂(t, q, p)) +
                      (db₁dx₂(t, q) * dv₂dq₃(t, q, p) + b₁(t, q) * d²v₂dq₂dq₃(t, q, p)) -
                      (db₂dx₂(t, q) * dv₁dq₃(t, q, p) + b₂(t, q) * d²v₁dq₂dq₃(t, q, p))

d²g₁dq₃dq₁(t, q, p) = (d²b₁dx₃dx₁(t, q) * v₂(t, q, p) + db₁dx₁(t, q) * dv₂dq₃(t, q, p)) -
                      (d²b₂dx₃dx₁(t, q) * v₁(t, q, p) + db₂dx₁(t, q) * dv₁dq₃(t, q, p)) +
                      (db₁dx₃(t, q) * dv₂dq₁(t, q, p) + b₁(t, q) * d²v₂dq₁dq₃(t, q, p)) -
                      (db₂dx₃(t, q) * dv₁dq₁(t, q, p) + b₂(t, q) * d²v₁dq₁dq₃(t, q, p))

d²g₁dq₃dq₂(t, q, p) = (d²b₁dx₃dx₂(t, q) * v₂(t, q, p) + db₁dx₂(t, q) * dv₂dq₃(t, q, p)) -
                      (d²b₂dx₃dx₂(t, q) * v₁(t, q, p) + db₂dx₂(t, q) * dv₁dq₃(t, q, p)) +
                      (db₁dx₃(t, q) * dv₂dq₂(t, q, p) + b₁(t, q) * d²v₂dq₂dq₃(t, q, p)) -
                      (db₂dx₃(t, q) * dv₁dq₂(t, q, p) + b₂(t, q) * d²v₁dq₂dq₃(t, q, p))

d²g₁dq₃dq₃(t, q, p) = (d²b₁dx₃dx₃(t, q) * v₂(t, q, p) + db₁dx₃(t, q) * dv₂dq₃(t, q, p)) -
                      (d²b₂dx₃dx₃(t, q) * v₁(t, q, p) + db₂dx₃(t, q) * dv₁dq₃(t, q, p)) +
                      (db₁dx₃(t, q) * dv₂dq₃(t, q, p) + b₁(t, q) * d²v₂dq₃dq₃(t, q, p)) -
                      (db₂dx₃(t, q) * dv₁dq₃(t, q, p) + b₂(t, q) * d²v₁dq₃dq₃(t, q, p))


d²g₂dq₁dq₁(t, q, p) = (d²b₁dx₁dx₁(t, q) * v₃(t, q, p) + db₁dx₁(t, q) * dv₃dq₁(t, q, p)) -
                      (d²b₃dx₁dx₁(t, q) * v₁(t, q, p) + db₃dx₁(t, q) * dv₁dq₁(t, q, p)) +
                      (db₁dx₁(t, q) * dv₃dq₁(t, q, p) + b₁(t, q) * d²v₃dq₁dq₁(t, q, p)) -
                      (db₃dx₁(t, q) * dv₁dq₁(t, q, p) + b₃(t, q) * d²v₁dq₁dq₁(t, q, p))

d²g₂dq₁dq₂(t, q, p) = (d²b₁dx₁dx₂(t, q) * v₃(t, q, p) + db₁dx₂(t, q) * dv₃dq₁(t, q, p)) -
                      (d²b₃dx₁dx₂(t, q) * v₁(t, q, p) + db₃dx₂(t, q) * dv₁dq₁(t, q, p)) +
                      (db₁dx₁(t, q) * dv₃dq₂(t, q, p) + b₁(t, q) * d²v₃dq₁dq₂(t, q, p)) -
                      (db₃dx₁(t, q) * dv₁dq₂(t, q, p) + b₃(t, q) * d²v₁dq₁dq₂(t, q, p))

d²g₂dq₁dq₃(t, q, p) = (d²b₁dx₁dx₃(t, q) * v₃(t, q, p) + db₁dx₃(t, q) * dv₃dq₁(t, q, p)) -
                      (d²b₃dx₁dx₃(t, q) * v₁(t, q, p) + db₃dx₃(t, q) * dv₁dq₁(t, q, p)) +
                      (db₁dx₁(t, q) * dv₃dq₃(t, q, p) + b₁(t, q) * d²v₃dq₁dq₃(t, q, p)) -
                      (db₃dx₁(t, q) * dv₁dq₃(t, q, p) + b₃(t, q) * d²v₁dq₁dq₃(t, q, p))

d²g₂dq₂dq₁(t, q, p) = (d²b₁dx₂dx₁(t, q) * v₃(t, q, p) + db₁dx₁(t, q) * dv₃dq₂(t, q, p)) -
                      (d²b₃dx₂dx₁(t, q) * v₁(t, q, p) + db₃dx₁(t, q) * dv₁dq₂(t, q, p)) +
                      (db₁dx₂(t, q) * dv₃dq₁(t, q, p) + b₁(t, q) * d²v₃dq₁dq₂(t, q, p)) -
                      (db₃dx₂(t, q) * dv₁dq₁(t, q, p) + b₃(t, q) * d²v₁dq₁dq₂(t, q, p))

d²g₂dq₂dq₂(t, q, p) = (d²b₁dx₂dx₂(t, q) * v₃(t, q, p) + db₁dx₂(t, q) * dv₃dq₂(t, q, p)) -
                      (d²b₃dx₂dx₂(t, q) * v₁(t, q, p) + db₃dx₂(t, q) * dv₁dq₂(t, q, p)) +
                      (db₁dx₂(t, q) * dv₃dq₂(t, q, p) + b₁(t, q) * d²v₃dq₂dq₂(t, q, p)) -
                      (db₃dx₂(t, q) * dv₁dq₂(t, q, p) + b₃(t, q) * d²v₁dq₂dq₂(t, q, p))

d²g₂dq₂dq₃(t, q, p) = (d²b₁dx₂dx₃(t, q) * v₃(t, q, p) + db₁dx₃(t, q) * dv₃dq₂(t, q, p)) -
                      (d²b₃dx₂dx₃(t, q) * v₁(t, q, p) + db₃dx₃(t, q) * dv₁dq₂(t, q, p)) +
                      (db₁dx₂(t, q) * dv₃dq₃(t, q, p) + b₁(t, q) * d²v₃dq₂dq₃(t, q, p)) -
                      (db₃dx₂(t, q) * dv₁dq₃(t, q, p) + b₃(t, q) * d²v₁dq₂dq₃(t, q, p))

d²g₂dq₃dq₁(t, q, p) = (d²b₁dx₃dx₁(t, q) * v₃(t, q, p) + db₁dx₁(t, q) * dv₃dq₃(t, q, p)) -
                      (d²b₃dx₃dx₁(t, q) * v₁(t, q, p) + db₃dx₁(t, q) * dv₁dq₃(t, q, p)) +
                      (db₁dx₃(t, q) * dv₃dq₁(t, q, p) + b₁(t, q) * d²v₃dq₁dq₃(t, q, p)) -
                      (db₃dx₃(t, q) * dv₁dq₁(t, q, p) + b₃(t, q) * d²v₁dq₁dq₃(t, q, p))

d²g₂dq₃dq₂(t, q, p) = (d²b₁dx₃dx₂(t, q) * v₃(t, q, p) + db₁dx₂(t, q) * dv₃dq₃(t, q, p)) -
                      (d²b₃dx₃dx₂(t, q) * v₁(t, q, p) + db₃dx₂(t, q) * dv₁dq₃(t, q, p)) +
                      (db₁dx₃(t, q) * dv₃dq₂(t, q, p) + b₁(t, q) * d²v₃dq₂dq₃(t, q, p)) -
                      (db₃dx₃(t, q) * dv₁dq₂(t, q, p) + b₃(t, q) * d²v₁dq₂dq₃(t, q, p))

d²g₂dq₃dq₃(t, q, p) = (d²b₁dx₃dx₃(t, q) * v₃(t, q, p) + db₁dx₃(t, q) * dv₃dq₃(t, q, p)) -
                      (d²b₃dx₃dx₃(t, q) * v₁(t, q, p) + db₃dx₃(t, q) * dv₁dq₃(t, q, p)) +
                      (db₁dx₃(t, q) * dv₃dq₃(t, q, p) + b₁(t, q) * d²v₃dq₃dq₃(t, q, p)) -
                      (db₃dx₃(t, q) * dv₁dq₃(t, q, p) + b₃(t, q) * d²v₁dq₃dq₃(t, q, p))


d²g₁dq₁dp₁(t, q, p) = db₁dx₁(t, q) * dv₂dp₁(t, q, p) - db₂dx₁(t, q) * dv₁dp₁(t, q, p)
d²g₁dq₁dp₂(t, q, p) = db₁dx₁(t, q) * dv₂dp₂(t, q, p) - db₂dx₁(t, q) * dv₁dp₂(t, q, p)
d²g₁dq₁dp₃(t, q, p) = db₁dx₁(t, q) * dv₂dp₃(t, q, p) - db₂dx₁(t, q) * dv₁dp₃(t, q, p)

d²g₁dq₂dp₁(t, q, p) = db₁dx₂(t, q) * dv₂dp₁(t, q, p) - db₂dx₂(t, q) * dv₁dp₁(t, q, p)
d²g₁dq₂dp₂(t, q, p) = db₁dx₂(t, q) * dv₂dp₂(t, q, p) - db₂dx₂(t, q) * dv₁dp₂(t, q, p)
d²g₁dq₂dp₃(t, q, p) = db₁dx₂(t, q) * dv₂dp₃(t, q, p) - db₂dx₂(t, q) * dv₁dp₃(t, q, p)

d²g₁dq₃dp₁(t, q, p) = db₁dx₃(t, q) * dv₂dp₁(t, q, p) - db₂dx₃(t, q) * dv₁dp₁(t, q, p)
d²g₁dq₃dp₂(t, q, p) = db₁dx₃(t, q) * dv₂dp₂(t, q, p) - db₂dx₃(t, q) * dv₁dp₂(t, q, p)
d²g₁dq₃dp₃(t, q, p) = db₁dx₃(t, q) * dv₂dp₃(t, q, p) - db₂dx₃(t, q) * dv₁dp₃(t, q, p)


d²g₂dq₁dp₁(t, q, p) = db₁dx₁(t, q) * dv₃dp₁(t, q, p) - db₃dx₁(t, q) * dv₁dp₁(t, q, p)
d²g₂dq₁dp₂(t, q, p) = db₁dx₁(t, q) * dv₃dp₂(t, q, p) - db₃dx₁(t, q) * dv₁dp₂(t, q, p)
d²g₂dq₁dp₃(t, q, p) = db₁dx₁(t, q) * dv₃dp₃(t, q, p) - db₃dx₁(t, q) * dv₁dp₃(t, q, p)

d²g₂dq₂dp₁(t, q, p) = db₁dx₂(t, q) * dv₃dp₁(t, q, p) - db₃dx₂(t, q) * dv₁dp₁(t, q, p)
d²g₂dq₂dp₂(t, q, p) = db₁dx₂(t, q) * dv₃dp₂(t, q, p) - db₃dx₂(t, q) * dv₁dp₂(t, q, p)
d²g₂dq₂dp₃(t, q, p) = db₁dx₂(t, q) * dv₃dp₃(t, q, p) - db₃dx₂(t, q) * dv₁dp₃(t, q, p)

d²g₂dq₃dp₁(t, q, p) = db₁dx₃(t, q) * dv₃dp₁(t, q, p) - db₃dx₃(t, q) * dv₁dp₁(t, q, p)
d²g₂dq₃dp₂(t, q, p) = db₁dx₃(t, q) * dv₃dp₂(t, q, p) - db₃dx₃(t, q) * dv₁dp₂(t, q, p)
d²g₂dq₃dp₃(t, q, p) = db₁dx₃(t, q) * dv₃dp₃(t, q, p) - db₃dx₃(t, q) * dv₁dp₃(t, q, p)


d²g₁dp₁dp₁(t, q, p) = zero(eltype(q))
d²g₁dp₁dp₂(t, q, p) = zero(eltype(q))
d²g₁dp₁dp₃(t, q, p) = zero(eltype(q))

d²g₁dp₂dp₁(t, q, p) = zero(eltype(q))
d²g₁dp₂dp₂(t, q, p) = zero(eltype(q))
d²g₁dp₂dp₃(t, q, p) = zero(eltype(q))

d²g₁dp₃dp₁(t, q, p) = zero(eltype(q))
d²g₁dp₃dp₂(t, q, p) = zero(eltype(q))
d²g₁dp₃dp₃(t, q, p) = zero(eltype(q))


d²g₂dp₁dp₁(t, q, p) = zero(eltype(q))
d²g₂dp₁dp₂(t, q, p) = zero(eltype(q))
d²g₂dp₁dp₃(t, q, p) = zero(eltype(q))

d²g₂dp₂dp₁(t, q, p) = zero(eltype(q))
d²g₂dp₂dp₂(t, q, p) = zero(eltype(q))
d²g₂dp₂dp₃(t, q, p) = zero(eltype(q))

d²g₂dp₃dp₁(t, q, p) = zero(eltype(q))
d²g₂dp₃dp₂(t, q, p) = zero(eltype(q))
d²g₂dp₃dp₃(t, q, p) = zero(eltype(q))



dλₒdq₁(t, q, p) = (
    d²g₁dq₁dq₁(t, q, p) * dg₂dp₁(t, q, p) +
    d²g₁dq₁dq₂(t, q, p) * dg₂dp₂(t, q, p) +
    d²g₁dq₁dq₃(t, q, p) * dg₂dp₃(t, q, p) -
    d²g₁dq₁dp₁(t, q, p) * dg₂dq₁(t, q, p) -
    d²g₁dq₁dp₂(t, q, p) * dg₂dq₂(t, q, p) -
    d²g₁dq₁dp₃(t, q, p) * dg₂dq₃(t, q, p) +
    dg₁dq₁(t, q, p) * d²g₂dq₁dp₁(t, q, p) +
    dg₁dq₂(t, q, p) * d²g₂dq₁dp₂(t, q, p) +
    dg₁dq₃(t, q, p) * d²g₂dq₁dp₃(t, q, p) -
    dg₁dp₁(t, q, p) * d²g₂dq₁dq₁(t, q, p) -
    dg₁dp₂(t, q, p) * d²g₂dq₁dq₂(t, q, p) -
    dg₁dp₃(t, q, p) * d²g₂dq₁dq₃(t, q, p)
)

dλₒdq₂(t, q, p) = (
    d²g₁dq₂dq₁(t, q, p) * dg₂dp₁(t, q, p) +
    d²g₁dq₂dq₂(t, q, p) * dg₂dp₂(t, q, p) +
    d²g₁dq₂dq₃(t, q, p) * dg₂dp₃(t, q, p) -
    d²g₁dq₂dp₁(t, q, p) * dg₂dq₁(t, q, p) -
    d²g₁dq₂dp₂(t, q, p) * dg₂dq₂(t, q, p) -
    d²g₁dq₂dp₃(t, q, p) * dg₂dq₃(t, q, p) +
    dg₁dq₁(t, q, p) * d²g₂dq₂dp₁(t, q, p) +
    dg₁dq₂(t, q, p) * d²g₂dq₂dp₂(t, q, p) +
    dg₁dq₃(t, q, p) * d²g₂dq₂dp₃(t, q, p) -
    dg₁dp₁(t, q, p) * d²g₂dq₁dq₂(t, q, p) -
    dg₁dp₂(t, q, p) * d²g₂dq₂dq₂(t, q, p) -
    dg₁dp₃(t, q, p) * d²g₂dq₃dq₂(t, q, p)
)

dλₒdq₃(t, q, p) = (
    d²g₁dq₃dq₁(t, q, p) * dg₂dp₁(t, q, p) +
    d²g₁dq₃dq₂(t, q, p) * dg₂dp₂(t, q, p) +
    d²g₁dq₃dq₃(t, q, p) * dg₂dp₃(t, q, p) -
    d²g₁dq₃dp₁(t, q, p) * dg₂dq₁(t, q, p) -
    d²g₁dq₃dp₂(t, q, p) * dg₂dq₂(t, q, p) -
    d²g₁dq₃dp₃(t, q, p) * dg₂dq₃(t, q, p) +
    dg₁dq₁(t, q, p) * d²g₂dq₃dp₁(t, q, p) +
    dg₁dq₂(t, q, p) * d²g₂dq₃dp₂(t, q, p) +
    dg₁dq₃(t, q, p) * d²g₂dq₃dp₃(t, q, p) -
    dg₁dp₁(t, q, p) * d²g₂dq₃dq₁(t, q, p) -
    dg₁dp₂(t, q, p) * d²g₂dq₃dq₂(t, q, p) -
    dg₁dp₃(t, q, p) * d²g₂dq₃dq₃(t, q, p)
)

dλₒdp₁(t, q, p) = (
    d²g₁dq₁dp₁(t, q, p) * dg₂dp₁(t, q, p) +
    d²g₁dq₂dp₁(t, q, p) * dg₂dp₂(t, q, p) +
    d²g₁dq₃dp₁(t, q, p) * dg₂dp₃(t, q, p) -
    d²g₁dp₁dp₁(t, q, p) * dg₂dq₁(t, q, p) -
    d²g₁dp₁dp₂(t, q, p) * dg₂dq₂(t, q, p) -
    d²g₁dp₁dp₃(t, q, p) * dg₂dq₃(t, q, p) +
    dg₁dq₁(t, q, p) * d²g₂dp₁dp₁(t, q, p) +
    dg₁dq₂(t, q, p) * d²g₂dp₂dp₁(t, q, p) +
    dg₁dq₃(t, q, p) * d²g₂dp₃dp₁(t, q, p) -
    dg₁dp₁(t, q, p) * d²g₂dq₁dp₁(t, q, p) -
    dg₁dp₂(t, q, p) * d²g₂dq₂dp₁(t, q, p) -
    dg₁dp₃(t, q, p) * d²g₂dq₃dp₁(t, q, p)
)

dλₒdp₂(t, q, p) = (
    d²g₁dq₁dp₂(t, q, p) * dg₂dp₁(t, q, p) +
    d²g₁dq₂dp₂(t, q, p) * dg₂dp₂(t, q, p) +
    d²g₁dq₃dp₂(t, q, p) * dg₂dp₃(t, q, p) -
    d²g₁dp₁dp₂(t, q, p) * dg₂dq₁(t, q, p) -
    d²g₁dp₂dp₂(t, q, p) * dg₂dq₂(t, q, p) -
    d²g₁dp₃dp₂(t, q, p) * dg₂dq₃(t, q, p) +
    dg₁dq₁(t, q, p) * d²g₂dp₁dp₂(t, q, p) +
    dg₁dq₂(t, q, p) * d²g₂dp₂dp₂(t, q, p) +
    dg₁dq₃(t, q, p) * d²g₂dp₃dp₂(t, q, p) -
    dg₁dp₁(t, q, p) * d²g₂dq₁dp₂(t, q, p) -
    dg₁dp₂(t, q, p) * d²g₂dq₂dp₂(t, q, p) -
    dg₁dp₃(t, q, p) * d²g₂dq₃dp₂(t, q, p)
)

dλₒdp₃(t, q, p) = (
    d²g₁dq₁dp₃(t, q, p) * dg₂dp₁(t, q, p) +
    d²g₁dq₂dp₃(t, q, p) * dg₂dp₂(t, q, p) +
    d²g₁dq₃dp₃(t, q, p) * dg₂dp₃(t, q, p) -
    d²g₁dp₁dp₃(t, q, p) * dg₂dq₁(t, q, p) -
    d²g₁dp₂dp₃(t, q, p) * dg₂dq₂(t, q, p) -
    d²g₁dp₃dp₃(t, q, p) * dg₂dq₃(t, q, p) +
    dg₁dq₁(t, q, p) * d²g₂dp₁dp₃(t, q, p) +
    dg₁dq₂(t, q, p) * d²g₂dp₂dp₃(t, q, p) +
    dg₁dq₃(t, q, p) * d²g₂dp₃dp₃(t, q, p) -
    dg₁dp₁(t, q, p) * d²g₂dq₁dp₃(t, q, p) -
    dg₁dp₂(t, q, p) * d²g₂dq₂dp₃(t, q, p) -
    dg₁dp₃(t, q, p) * d²g₂dq₃dp₃(t, q, p)
)

dλ₁dq₁(t, q, p, params) = -dλₒdq₁(t, q, p) * λ₁(t, q, p, params) / λₒ(t, q, p) + (
    d²g₂dq₁dq₁(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₂dq₁dq₂(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₂dq₁dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₂dq₁dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₂dq₁dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₂dq₁dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₂dq₁(t, q, p) * d²Hdq₁dp₁(t, q, p, params) +
    dg₂dq₂(t, q, p) * d²Hdq₁dp₂(t, q, p, params) +
    dg₂dq₃(t, q, p) * d²Hdq₁dp₃(t, q, p, params) -
    dg₂dp₁(t, q, p) * d²Hdq₁dq₁(t, q, p, params) -
    dg₂dp₂(t, q, p) * d²Hdq₁dq₂(t, q, p, params) -
    dg₂dp₃(t, q, p) * d²Hdq₁dq₃(t, q, p, params)
) / λₒ(t, q, p)

dλ₁dq₂(t, q, p, params) = -dλₒdq₂(t, q, p) * λ₁(t, q, p, params) / λₒ(t, q, p) + (
    d²g₂dq₁dq₂(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₂dq₂dq₂(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₂dq₂dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₂dq₂dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₂dq₂dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₂dq₂dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₂dq₁(t, q, p) * d²Hdq₂dp₁(t, q, p, params) +
    dg₂dq₂(t, q, p) * d²Hdq₂dp₂(t, q, p, params) +
    dg₂dq₃(t, q, p) * d²Hdq₂dp₃(t, q, p, params) -
    dg₂dp₁(t, q, p) * d²Hdq₁dq₂(t, q, p, params) -
    dg₂dp₂(t, q, p) * d²Hdq₂dq₂(t, q, p, params) -
    dg₂dp₃(t, q, p) * d²Hdq₂dq₃(t, q, p, params)
) / λₒ(t, q, p)

dλ₁dq₃(t, q, p, params) = -dλₒdq₃(t, q, p) * λ₁(t, q, p, params) / λₒ(t, q, p) + (
    d²g₂dq₁dq₃(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₂dq₂dq₃(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₂dq₃dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₂dq₃dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₂dq₃dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₂dq₃dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₂dq₁(t, q, p) * d²Hdq₃dp₁(t, q, p, params) +
    dg₂dq₂(t, q, p) * d²Hdq₃dp₂(t, q, p, params) +
    dg₂dq₃(t, q, p) * d²Hdq₃dp₃(t, q, p, params) -
    dg₂dp₁(t, q, p) * d²Hdq₁dq₃(t, q, p, params) -
    dg₂dp₂(t, q, p) * d²Hdq₂dq₃(t, q, p, params) -
    dg₂dp₃(t, q, p) * d²Hdq₃dq₃(t, q, p, params)
) / λₒ(t, q, p)

dλ₁dp₁(t, q, p, params) = -dλₒdp₁(t, q, p) * λ₁(t, q, p, params) / λₒ(t, q, p) + (
    d²g₂dq₁dp₁(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₂dq₂dp₁(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₂dq₃dp₁(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₂dp₁dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₂dp₁dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₂dp₁dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₂dq₁(t, q, p) * d²Hdp₁dp₁(t, q, p, params) +
    dg₂dq₂(t, q, p) * d²Hdp₁dp₂(t, q, p, params) +
    dg₂dq₃(t, q, p) * d²Hdp₁dp₃(t, q, p, params) -
    dg₂dp₁(t, q, p) * d²Hdq₁dp₁(t, q, p, params) -
    dg₂dp₂(t, q, p) * d²Hdq₂dp₁(t, q, p, params) -
    dg₂dp₃(t, q, p) * d²Hdq₃dp₁(t, q, p, params)
) / λₒ(t, q, p)

dλ₁dp₂(t, q, p, params) = -dλₒdp₂(t, q, p) * λ₁(t, q, p, params) / λₒ(t, q, p) + (
    d²g₂dq₁dp₂(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₂dq₂dp₂(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₂dq₃dp₂(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₂dp₁dp₂(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₂dp₂dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₂dp₂dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₂dq₁(t, q, p) * d²Hdp₁dp₂(t, q, p, params) +
    dg₂dq₂(t, q, p) * d²Hdp₂dp₂(t, q, p, params) +
    dg₂dq₃(t, q, p) * d²Hdp₂dp₃(t, q, p, params) -
    dg₂dp₁(t, q, p) * d²Hdq₁dp₂(t, q, p, params) -
    dg₂dp₂(t, q, p) * d²Hdq₂dp₂(t, q, p, params) -
    dg₂dp₃(t, q, p) * d²Hdq₃dp₂(t, q, p, params)
) / λₒ(t, q, p)

dλ₁dp₃(t, q, p, params) = -dλₒdp₃(t, q, p) * λ₁(t, q, p, params) / λₒ(t, q, p) + (
    d²g₂dq₁dp₃(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₂dq₂dp₃(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₂dq₃dp₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₂dp₁dp₃(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₂dp₂dp₃(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₂dp₃dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₂dq₁(t, q, p) * d²Hdp₁dp₃(t, q, p, params) +
    dg₂dq₂(t, q, p) * d²Hdp₂dp₃(t, q, p, params) +
    dg₂dq₃(t, q, p) * d²Hdp₃dp₃(t, q, p, params) -
    dg₂dp₁(t, q, p) * d²Hdq₁dp₃(t, q, p, params) -
    dg₂dp₂(t, q, p) * d²Hdq₂dp₃(t, q, p, params) -
    dg₂dp₃(t, q, p) * d²Hdq₃dp₃(t, q, p, params)
) / λₒ(t, q, p)

dλ₂dq₁(t, q, p, params) = +dλₒdq₁(t, q, p) * λ₂(t, q, p, params) / λₒ(t, q, p) - (
    d²g₁dq₁dq₁(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₁dq₁dq₂(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₁dq₁dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₁dq₁dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₁dq₁dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₁dq₁dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₁dq₁(t, q, p) * d²Hdq₁dp₁(t, q, p, params) +
    dg₁dq₂(t, q, p) * d²Hdq₁dp₂(t, q, p, params) +
    dg₁dq₃(t, q, p) * d²Hdq₁dp₃(t, q, p, params) -
    dg₁dp₁(t, q, p) * d²Hdq₁dq₁(t, q, p, params) -
    dg₁dp₂(t, q, p) * d²Hdq₁dq₂(t, q, p, params) -
    dg₁dp₃(t, q, p) * d²Hdq₁dq₃(t, q, p, params)
) / λₒ(t, q, p)

dλ₂dq₂(t, q, p, params) = +dλₒdq₂(t, q, p) * λ₂(t, q, p, params) / λₒ(t, q, p) - (
    d²g₁dq₁dq₂(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₁dq₂dq₂(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₁dq₂dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₁dq₂dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₁dq₂dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₁dq₂dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₁dq₁(t, q, p) * d²Hdq₂dp₁(t, q, p, params) +
    dg₁dq₂(t, q, p) * d²Hdq₂dp₂(t, q, p, params) +
    dg₁dq₃(t, q, p) * d²Hdq₂dp₃(t, q, p, params) -
    dg₁dp₁(t, q, p) * d²Hdq₁dq₂(t, q, p, params) -
    dg₁dp₂(t, q, p) * d²Hdq₂dq₂(t, q, p, params) -
    dg₁dp₃(t, q, p) * d²Hdq₂dq₃(t, q, p, params)
) / λₒ(t, q, p)

dλ₂dq₃(t, q, p, params) = +dλₒdq₃(t, q, p) * λ₂(t, q, p, params) / λₒ(t, q, p) - (
    d²g₁dq₁dq₃(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₁dq₂dq₃(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₁dq₃dq₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₁dq₃dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₁dq₃dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₁dq₃dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₁dq₁(t, q, p) * d²Hdq₃dp₁(t, q, p, params) +
    dg₁dq₂(t, q, p) * d²Hdq₃dp₂(t, q, p, params) +
    dg₁dq₃(t, q, p) * d²Hdq₃dp₃(t, q, p, params) -
    dg₁dp₁(t, q, p) * d²Hdq₁dq₃(t, q, p, params) -
    dg₁dp₂(t, q, p) * d²Hdq₂dq₃(t, q, p, params) -
    dg₁dp₃(t, q, p) * d²Hdq₃dq₃(t, q, p, params)
) / λₒ(t, q, p)

dλ₂dp₁(t, q, p, params) = +dλₒdp₁(t, q, p) * λ₂(t, q, p, params) / λₒ(t, q, p) - (
    d²g₁dq₁dp₁(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₁dq₂dp₁(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₁dq₃dp₁(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₁dp₁dp₁(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₁dp₁dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₁dp₁dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₁dq₁(t, q, p) * d²Hdp₁dp₁(t, q, p, params) +
    dg₁dq₂(t, q, p) * d²Hdp₁dp₂(t, q, p, params) +
    dg₁dq₃(t, q, p) * d²Hdp₁dp₃(t, q, p, params) -
    dg₁dp₁(t, q, p) * d²Hdq₁dp₁(t, q, p, params) -
    dg₁dp₂(t, q, p) * d²Hdq₂dp₁(t, q, p, params) -
    dg₁dp₃(t, q, p) * d²Hdq₃dp₁(t, q, p, params)
) / λₒ(t, q, p)

dλ₂dp₂(t, q, p, params) = +dλₒdp₂(t, q, p) * λ₂(t, q, p, params) / λₒ(t, q, p) - (
    d²g₁dq₁dp₂(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₁dq₂dp₂(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₁dq₃dp₂(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₁dp₁dp₂(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₁dp₂dp₂(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₁dp₂dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₁dq₁(t, q, p) * d²Hdp₁dp₂(t, q, p, params) +
    dg₁dq₂(t, q, p) * d²Hdp₂dp₂(t, q, p, params) +
    dg₁dq₃(t, q, p) * d²Hdp₂dp₃(t, q, p, params) -
    dg₁dp₁(t, q, p) * d²Hdq₁dp₂(t, q, p, params) -
    dg₁dp₂(t, q, p) * d²Hdq₂dp₂(t, q, p, params) -
    dg₁dp₃(t, q, p) * d²Hdq₃dp₂(t, q, p, params)
) / λₒ(t, q, p)

dλ₂dp₃(t, q, p, params) = +dλₒdp₃(t, q, p) * λ₂(t, q, p, params) / λₒ(t, q, p) - (
    d²g₁dq₁dp₃(t, q, p) * dHdp₁(t, q, p, params) +
    d²g₁dq₂dp₃(t, q, p) * dHdp₂(t, q, p, params) +
    d²g₁dq₃dp₃(t, q, p) * dHdp₃(t, q, p, params) -
    d²g₁dp₁dp₃(t, q, p) * dHdq₁(t, q, p, params) -
    d²g₁dp₂dp₃(t, q, p) * dHdq₂(t, q, p, params) -
    d²g₁dp₃dp₃(t, q, p) * dHdq₃(t, q, p, params) +
    dg₁dq₁(t, q, p) * d²Hdp₁dp₃(t, q, p, params) +
    dg₁dq₂(t, q, p) * d²Hdp₂dp₃(t, q, p, params) +
    dg₁dq₃(t, q, p) * d²Hdp₃dp₃(t, q, p, params) -
    dg₁dp₁(t, q, p) * d²Hdq₁dp₃(t, q, p, params) -
    dg₁dp₂(t, q, p) * d²Hdq₂dp₃(t, q, p, params) -
    dg₁dp₃(t, q, p) * d²Hdq₃dp₃(t, q, p, params)
) / λₒ(t, q, p)



function guiding_center_3d_canonical_v(v, t, q, p, params)
    guiding_center_3d_v(v, t, q, p, params)

    v[1] += dλ₁dp₁(t, q, p, params) * g₁(t, q, p) +
            dλ₂dp₁(t, q, p, params) * g₂(t, q, p)
    v[2] += dλ₁dp₂(t, q, p, params) * g₁(t, q, p) +
            dλ₂dp₂(t, q, p, params) * g₂(t, q, p)
    v[3] += dλ₁dp₃(t, q, p, params) * g₁(t, q, p) +
            dλ₂dp₃(t, q, p, params) * g₂(t, q, p)

    nothing
end

function guiding_center_3d_canonical_f(f, t, q, p, params)
    guiding_center_3d_f(f, t, q, p, params)

    f[1] -= dλ₁dq₁(t, q, p, params) * g₁(t, q, p) +
            dλ₂dq₁(t, q, p, params) * g₂(t, q, p)
    f[2] -= dλ₁dq₂(t, q, p, params) * g₁(t, q, p) +
            dλ₂dq₂(t, q, p, params) * g₂(t, q, p)
    f[3] -= dλ₁dq₃(t, q, p, params) * g₁(t, q, p) +
            dλ₂dq₃(t, q, p, params) * g₂(t, q, p)

    nothing
end


function hode_canonical(q₀, p₀, parameters; tspan=tspan, tstep=Δt, periodic=true)
    # println("3D Guiding Center model initial constraints g₁ = $(g₁(t₀, q₀, p₀)) and g₂ = $(g₂(t₀, q₀, p₀))")
    HODEProblem(
        guiding_center_3d_canonical_v,
        guiding_center_3d_canonical_f,
        hamiltonian_canonical,
        tspan, tstep, q₀, p₀;
        parameters=parameters,
        periodicity=guiding_center_3d_periodicity(q₀, periodic))
end

function hode_canonical(x₀, parameters; tspan=tspan, kwargs...)
    hode_canonical(initial_conditions(tspan[begin], x₀)..., parameters; tspan=tspan, kwargs...)
end
