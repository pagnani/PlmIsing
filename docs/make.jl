using Documenter, PlmIsing

push!(LOAD_PATH,"../src/")
makedocs(sitename = "PlmIsing",
	 modules = [PlmIsing],
	 doctest = true)
deploydocs(
	   branch = "gh-pages",
	   repo = "github.com/pagnani/PlmIsing.git",
	   versions = ["stable" => "v^"]
	  )
