workspace "Ionizer"
	architecture "x64"
	configurations{
		"Debug",
		"Release",
		"Dist"
	}

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"

project "Ionizer"
	location "Ionizer"
	kind "ConsoleApp"
	language "C++"
	
	openmp "On"
	
	targetdir ("bin/" .. outputdir .. "/%{prj.name}")
	objdir ("bin-int/" .. outputdir .. "/%{prj.name}")
	
	files {
		"%{prj.name}/src/**.h",
		"%{prj.name}/src/**.cpp"
	}
	
	includedirs {
		"%{prj.name}/vendor/spdlog/include",
		"%{prj.name}/vendor/eigen"
	}
	
	filter "system:windows"
		cppdialect "C++17"
		staticruntime "On"
		systemversion "latest"
		
	filter "configurations:Debug"
		defines "ION_DEBUG"
		symbols "On"
		
	filter "configurations:Release"
		defines "ION_RELEASE"
		optimize "On"
		
	filter "configurations:Dist"
		defines "ION_DIST"
		optimize "On"
		
project "Sandbox"
	location "Sandbox"
	kind "ConsoleApp"
	language "C++"
	
	openmp "On"
	
	targetdir ("bin/" .. outputdir .. "/%{prj.name}")
	objdir ("bin-int/" .. outputdir .. "/%{prj.name}")
	
	files {
		"%{prj.name}/src/**.h",
		"%{prj.name}/src/**.cpp"
	}
	
	includedirs {
		"Ionizer/src"
	}
	
	filter "system:windows"
		cppdialect "C++17"
		staticruntime "On"
		systemversion "latest"
		
		defines {
			"ION_PLATFORM_WINDOWS"
		}
		
	filter "configurations:Debug"
		defines "ION_DEBUG"
		symbols "On"
		
	filter "configurations:Release"
		defines "ION_RELEASE"
		optimize "On"
		
	filter "configurations:Dist"
		defines "ION_DIST"
		optimize "On"