project "Ionizer"
	kind "StaticLib"
	language "C++"
	cppdialect "C++17"
	openmp "On"
	
	targetdir ("%{wks.location}/bin/" .. outputdir .. "/%{prj.name}")
	objdir ("%{wks.location}/bin-int/" .. outputdir .. "/%{prj.name}")
	
	pchheader "Ionpch.h"
	pchsource "src/Ionpch.cpp"
	
	files {
		"src/**.h",
		"src/**.cpp"
	}
	
	includedirs {
		"src",
		"vendor/spdlog/include",
		"vendor/eigen"
	}
	
	filter "system:windows"
		cppdialect "C++17"
		systemversion "latest"
		defines {
			"IONIZER_PLATFORM_WINDOWS",
			"IONIZER_BUILD_DLL"
		}
		
	filter "configurations:Debug"
		defines "IONIZER_DEBUG"
		runtime "Debug"
		symbols "On"
		
	filter "configurations:Release"
		defines "IONIZER_RELEASE"
		runtime "Release"
		optimize "On"
		
	filter "configurations:Dist"
		defines "IONIZER_DIST"
		runtime "Release"
		optimize "On"
	