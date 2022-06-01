workspace "Ionizer"
	architecture "x64"
	configurations{
		"Debug",
		"Release",
		"Dist"
	}
	startproject "Sandbox"

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"

include "Ionizer"
include "Sandbox"