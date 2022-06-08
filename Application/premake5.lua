include "vendor/Walnut/WalnutExternal.lua"

project "Application"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++17"
	
	files {
		"src/**.h",
		"src/**.cpp"
	}
	
	includedirs 
	{
		"%{wks.location}/Ionizer/src",
		"%{wks.location}/Ionizer/vendor/spdlog/include",
		"%{wks.location}/Ionizer/vendor/eigen",
		"vendor/Walnut/vendor/imgui",
		"vendor/Walnut/vendor/glfw/include",
		"vendor/Walnut/vendor/glm",
		"vendor/Walnut/Walnut/src",
		"%{IncludeDir.VulkanSDK}"
	}
	
	links
	{
	   "Walnut",
	   "Ionizer"
	}
	
	targetdir ("%{wks.location}/bin/" .. outputdir .. "/%{prj.name}")
	objdir ("%{wks.location}/bin-int/" .. outputdir .. "/%{prj.name}")
	
	filter "system:windows"
		systemversion "latest"
		defines { 
			"WL_PLATFORM_WINDOWS",
			"IONIZER_PLATFORM_WINDOWS"
		}

	filter "configurations:Debug"
		defines { 
			"WL_DEBUG",
			"IONIZER_DEBUG"
		}
		runtime "Debug"
		symbols "On"

	filter "configurations:Release"
		defines { 
			"WL_RELEASE",
			"IONIZER_RELEASE"
		}
		runtime "Release"
		optimize "On"
		symbols "On"

	filter "configurations:Dist"
		kind "WindowedApp"
		defines { 
			"WL_DIST",
			"IONIZER_DIST"
		}
		runtime "Release"
		optimize "On"
		symbols "Off"