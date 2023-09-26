### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f4ac78d5-c204-43df-a140-c0d51fd7ebaa
begin
	using Colors, ColorVectorSpace, ImageShow, FileIO, ImageIO
	using CommonMark
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ 446d5f05-0643-4e70-9a0b-5fb2f85e2fce
md"## 0. Setup"

# ╔═╡ 6cd4456c-5c1b-11ee-3a8b-eb49a1fb1869
md"## 1. Most of what we do in numerical computing is just taking advantage of calculus."

# ╔═╡ 7b446c3e-411f-4f91-b202-3385e004c5cd
md"### Derivatives"

# ╔═╡ 779fe839-130b-4df4-83e4-7bd903569662
cm"""

Recall that

<div align="center">

``\dfrac{d}{dx} f(x) = \lim\limits_{\Delta x \rightarrow 0}\dfrac{f(x + \Delta x) - f(x)}{\Delta x}``.

</div>

"""

# ╔═╡ cdba1db8-d5d0-4e6a-a75f-d09c62a895f6
md"""
Let's see this in action! Choose a function:
"""

# ╔═╡ 99c2ffa0-d2c5-4aa7-8fdc-4112f4b955b4
md"We know that"

# ╔═╡ 85a182ed-ad87-4abd-ab52-b74d565d7b52
md"so we can get"

# ╔═╡ 526899b1-4e20-4e3f-b9c4-2c5425a21eb9
md"""Now drag the slider below to change the value of ``\text{exp}`` and ``\Delta x = 10^{\text{exp}}``:"""

# ╔═╡ f5d4935a-bf64-4938-809d-2b74fa2e1541
@bind exp Slider(-7:0, show_value=true)

# ╔═╡ 471a55ab-69b8-4232-acc9-0f3f895dcad0
L"""

\text{exp} = %$(exp),~\Delta x = 10^{%$(exp)}

"""

# ╔═╡ 40c1b434-8f2c-4c84-a573-3b8941ec6627
md"We then have"

# ╔═╡ 0a73d6ae-3be4-47ab-b1bc-94259bce7e97
md"### Integrals"

# ╔═╡ e66561fb-10d5-41e0-aba5-20811301831c
cm"""

Recall that

<div align="center">

``\displaystyle\int_{a}^{b} f(x)~dx = \lim\limits_{n \rightarrow \infty}\sum\limits_{i = 1}^{n} f(x_{i})\Delta x``,

</div>

"""

# ╔═╡ c626afe6-57e6-4625-9ac7-e7f9cb5daae3
md"where"

# ╔═╡ 0c40ac6e-dbbd-4620-9cab-d0c5c43991f0
L"""
\Delta x = \dfrac{b - a}{n}
"""

# ╔═╡ 48cc5e4d-d40f-4a02-9263-4ce6373fbc4d
md"and"

# ╔═╡ 5b3611e3-cab1-46d1-96b2-3d334f5852b2
L"""
x_{i} = a + \left(\dfrac{b-a}{n}\right)i.
"""

# ╔═╡ de2685b8-faec-4993-9739-2d3e9f335436
md"""
Let's see this in action! Choose a function:
"""

# ╔═╡ 42ceee55-c13d-46ee-bfca-c679452753d6
md"We know that"

# ╔═╡ 0f828f1b-d374-4ce2-bc38-496e7725fed6
md"so we can get"

# ╔═╡ a2668bfa-2dd2-46f4-8431-950b34819fcf
md"""Now drag the slider below to change the value of ``\text{exp}`` and ``n = 10^{\text{exp}}``:"""

# ╔═╡ cb94ddfa-172e-4185-9d51-2d4823dec0f2
@bind exp₂ Slider(0:7, show_value=true)

# ╔═╡ 5378c49e-d79c-4109-9f39-8832deab2b3e
L"""

\text{exp} = %$(exp₂),~n = 10^{%$(exp₂)}

"""

# ╔═╡ c9a831f7-fb86-4c94-b724-068880d75d09
md"We then have"

# ╔═╡ 4a717a47-c1ed-48d4-a87b-2f449975470d
md"### LinAlg example"

# ╔═╡ c95fbc3e-a11a-4e4a-8305-6eda9f437df8
md"## 2. Matrices are linear transformations."

# ╔═╡ 4d55e733-a8dd-444a-92b8-cfca6f814cec
md"""
The elements of the ``2 \times 2`` matrix below are *scrubbable*; watch the values of ``L`` and ``\text{det}`` change as you edit the individual values!
"""

# ╔═╡ 34abae52-63ee-468f-89fd-1e18d1daa25d
begin
	scrubbable_range = -5:0.1:5

	md"""
	``(``	
	 $(@bind a Scrubbable(scrubbable_range; default=1.0))
	 $(@bind b Scrubbable(scrubbable_range; default=0.0))
	``)``
	
	``(``
	$(@bind c Scrubbable(scrubbable_range; default=0.0))
	$(@bind d Scrubbable(scrubbable_range; default=1.0))
	``)``
	"""
end

# ╔═╡ f933fabc-a8cf-415e-a52c-fa7b91708fb2
L"""

L = \begin{bmatrix}
%$(a) & %$(b) \\
%$(c) & %$(d)
\end{bmatrix}

,~\text{det} = %$(round(a*d - b*c; digits = 4))

"""

# ╔═╡ 68e245fd-7014-4545-9d9f-13ee791e7435
md"""
## Appendix

Behind-the-scenes stuff!
"""

# ╔═╡ d160ece5-93e4-4319-b2ae-c7b21129eccd
transform(a, b, c, d) = ((x₁, x₂),) -> [a*x₁ + b*x₂; c*x₁ + d*x₂]

# ╔═╡ 6f38c764-d3f8-4a51-8a4a-3b7702638185
derivative_sources = [
	(x -> 1, x -> 0, "1", "0") => "f(x) = 1",
	(x -> x, x -> 1, "x", "1") => "f(x) = x",
	(x -> x^2, x -> 2x, "x^{2}", "2x") => "f(x) = x²",
	(x -> log(x), x -> 1/x, "\\ln(x)", "\\dfrac{1}{x}") => "f(x) = ln(x)",
	(x -> sin(x), x -> cos(x), "\\sin(x)", "\\cos(x)") => "f(x) = sin(x)",
	(x -> cos(x), x -> -sin(x), "\\cos(x)", "-\\sin(x)") => "f(x) = cos(x)"
]

# ╔═╡ 7def1388-9c49-4bb4-bd27-5c758f0db02e
md"""
$(@bind derivative_selection Select(derivative_sources))
"""

# ╔═╡ 1de1ccbf-06d3-424a-8e7a-19216cb241ee
L"""
f(x) = %$(derivative_selection[3])
"""

# ╔═╡ 0e62bfdb-3041-45ac-a144-9dfa43beec28
L"""
f'(x) = %$(derivative_selection[4]),
"""

# ╔═╡ e86f80a5-75c3-4139-8d8e-7a8ef0296c83
L"""
f'(1) \approx %$(round(derivative_selection[2](1); digits = 6))~\text{(analytically).}
"""

# ╔═╡ aeee7d88-eb22-4adf-a611-48e44fcc8d55
integral_sources = [
	(x -> 1, x -> x, "1", "x") => "f(x) = 1",
	(x -> x, x -> x^2/2, "x", "\\dfrac{x^{2}}{2}") => "f(x) = x",
	(x -> x^2, x -> x^3/3, "x^{2}", "\\dfrac{x^{3}}{3}") => "f(x) = x²",
	(x -> 1/x, x -> log(x), "\\dfrac{1}{x}", "\\ln(x)") => "f(x) = 1/x",
	(x -> sin(x), x -> -cos(x), "\\sin(x)", "-\\cos(x)") => "f(x) = sin(x)",
	(x -> cos(x), x -> sin(x), "\\cos(x)", "\\sin(x)") => "f(x) = cos(x)"
]

# ╔═╡ 3fa1143c-f8ac-4ae8-8541-4ee9adba5322
md"""
$(@bind integral_selection Select(integral_sources))
"""

# ╔═╡ dc4cea1e-4920-4336-bc4a-a43532815f1f
L"""
f(x) = %$(integral_selection[3])
"""

# ╔═╡ 9ddfb070-70ec-4e2a-9dae-c981622f4fc2
L"""
\displaystyle\int f(x)~dx = %$(integral_selection[4]),
"""

# ╔═╡ 2ca18fed-5c42-49f3-8f60-60095b6077d2
L"""
\displaystyle\int_{1}^{2} f(x)~dx \approx %$(round(integral_selection[2](2) - integral_selection[2](1); digits = 6))~\text{(analytically).}
"""

# ╔═╡ 71a8b452-b192-475c-a270-b6b77e2dea20
Δx = 10.0^exp

# ╔═╡ 8bd4d340-c8db-4dc3-9f7e-40b69060e8fd
L"""
f'(1) \approx \dfrac{f(1 + 10^{%$(exp)}) - f(1)}{10^{%$(exp)}} \approx %$(round((derivative_selection[1](1+Δx) - derivative_selection[1](1))/Δx; digits = 6))~\text{(numerically).}
"""

# ╔═╡ 063fa03d-1b70-480c-93e2-48fc275d0da7
n = 10^exp₂

# ╔═╡ f276b4fa-1226-4a57-8c88-c25c29dcbf8f
L"""
\displaystyle\int_{1}^{2}f(x)~dx \approx \sum\limits_{i=1}^{10^{%$(exp₂)}}f\left(1 + \dfrac{i}{n}\right)~\dfrac{1}{n} \approx %$(round(sum([integral_selection[1](1 + i/n) for i in 1:n])/n; digits = 6))~\text{(numerically).}
"""

# ╔═╡ d3ef7ff3-ae0b-4967-bf32-390b46af4f1b
img_sources = [
	"https://github.com/daryll-ko/cs138/blob/main/assets/cartesian_plane.png?raw=true" => "Cartesian Plane",
	"https://github.com/daryll-ko/cs138/blob/main/assets/dcs_logo.png?raw=true" => "DCS Logo",
	"https://github.com/daryll-ko/cs138/blob/main/assets/scl_logo.png?raw=true" => "SCL Logo"
]

# ╔═╡ 0dcea395-8b61-4dfe-b7c5-a4c17380d766
md"""
Choose an image to play with!

$(@bind img_source Select(img_sources))
"""

# ╔═╡ 484d83e0-423a-449a-b982-04bcbcc33a66
img_source

# ╔═╡ 12e19726-8cae-4db5-b1bd-daec9a9c21bd
img_filename = download(img_source)

# ╔═╡ 6bffce36-e6ba-472c-a9e5-e46ae473e2df
vanilla_img = load(img_filename);

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ColorVectorSpace = "c3611d14-8923-5661-9e6a-0046d554d3a4"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
ImageIO = "82e4d734-157c-48bb-816b-45c225c6df19"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
ColorVectorSpace = "~0.10.0"
Colors = "~0.12.10"
CommonMark = "~0.8.12"
FileIO = "~1.16.1"
ImageIO = "~0.6.7"
ImageShow = "~0.3.8"
LaTeXStrings = "~1.3.0"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "6b530cd68ce25614939c9c3a5562212cf3e700dd"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "fc5d1d3443a124fde6e92d0260cd9e064eba69f8"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.1"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "8e59ea773deee525c99a8018409f64f19fb719e6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.7"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "327713faef2a3e5c80f96bf38d1fa26f7a6ae29e"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "9b02b27ac477cad98114584ff964e3052f656a0f"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "b7dc44cb005a7ef743b8fe98970afef003efdce7"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.6"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─446d5f05-0643-4e70-9a0b-5fb2f85e2fce
# ╠═f4ac78d5-c204-43df-a140-c0d51fd7ebaa
# ╟─6cd4456c-5c1b-11ee-3a8b-eb49a1fb1869
# ╟─7b446c3e-411f-4f91-b202-3385e004c5cd
# ╟─779fe839-130b-4df4-83e4-7bd903569662
# ╟─cdba1db8-d5d0-4e6a-a75f-d09c62a895f6
# ╟─7def1388-9c49-4bb4-bd27-5c758f0db02e
# ╟─1de1ccbf-06d3-424a-8e7a-19216cb241ee
# ╟─99c2ffa0-d2c5-4aa7-8fdc-4112f4b955b4
# ╟─0e62bfdb-3041-45ac-a144-9dfa43beec28
# ╟─85a182ed-ad87-4abd-ab52-b74d565d7b52
# ╟─e86f80a5-75c3-4139-8d8e-7a8ef0296c83
# ╟─526899b1-4e20-4e3f-b9c4-2c5425a21eb9
# ╟─f5d4935a-bf64-4938-809d-2b74fa2e1541
# ╟─471a55ab-69b8-4232-acc9-0f3f895dcad0
# ╟─40c1b434-8f2c-4c84-a573-3b8941ec6627
# ╟─8bd4d340-c8db-4dc3-9f7e-40b69060e8fd
# ╟─0a73d6ae-3be4-47ab-b1bc-94259bce7e97
# ╟─e66561fb-10d5-41e0-aba5-20811301831c
# ╟─c626afe6-57e6-4625-9ac7-e7f9cb5daae3
# ╟─0c40ac6e-dbbd-4620-9cab-d0c5c43991f0
# ╟─48cc5e4d-d40f-4a02-9263-4ce6373fbc4d
# ╟─5b3611e3-cab1-46d1-96b2-3d334f5852b2
# ╟─de2685b8-faec-4993-9739-2d3e9f335436
# ╟─3fa1143c-f8ac-4ae8-8541-4ee9adba5322
# ╟─dc4cea1e-4920-4336-bc4a-a43532815f1f
# ╟─42ceee55-c13d-46ee-bfca-c679452753d6
# ╟─9ddfb070-70ec-4e2a-9dae-c981622f4fc2
# ╟─0f828f1b-d374-4ce2-bc38-496e7725fed6
# ╟─2ca18fed-5c42-49f3-8f60-60095b6077d2
# ╟─a2668bfa-2dd2-46f4-8431-950b34819fcf
# ╟─cb94ddfa-172e-4185-9d51-2d4823dec0f2
# ╟─5378c49e-d79c-4109-9f39-8832deab2b3e
# ╟─c9a831f7-fb86-4c94-b724-068880d75d09
# ╟─f276b4fa-1226-4a57-8c88-c25c29dcbf8f
# ╟─4a717a47-c1ed-48d4-a87b-2f449975470d
# ╟─c95fbc3e-a11a-4e4a-8305-6eda9f437df8
# ╟─4d55e733-a8dd-444a-92b8-cfca6f814cec
# ╟─34abae52-63ee-468f-89fd-1e18d1daa25d
# ╟─f933fabc-a8cf-415e-a52c-fa7b91708fb2
# ╟─0dcea395-8b61-4dfe-b7c5-a4c17380d766
# ╟─68e245fd-7014-4545-9d9f-13ee791e7435
# ╟─d160ece5-93e4-4319-b2ae-c7b21129eccd
# ╟─6f38c764-d3f8-4a51-8a4a-3b7702638185
# ╟─aeee7d88-eb22-4adf-a611-48e44fcc8d55
# ╟─71a8b452-b192-475c-a270-b6b77e2dea20
# ╟─063fa03d-1b70-480c-93e2-48fc275d0da7
# ╟─d3ef7ff3-ae0b-4967-bf32-390b46af4f1b
# ╟─484d83e0-423a-449a-b982-04bcbcc33a66
# ╟─12e19726-8cae-4db5-b1bd-daec9a9c21bd
# ╟─6bffce36-e6ba-472c-a9e5-e46ae473e2df
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
