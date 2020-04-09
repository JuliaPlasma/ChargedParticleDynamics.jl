
using Weave


weave("guiding_center_dynamics_splitting.jmd",
  out_path="build",
  doctype = "github")

weave("guiding_center_dynamics_vprk.jmd",
  out_path="build",
  doctype = "github")

