#ifndef _GLUTILS_H_
#define _GLUTILS_H_

#include <GL/gl.h>
#include <GL/glu.h>
#include <map>
#include <string>

namespace RenderEngine
{
	class GLRenderSetup
	{
		typedef void ( *queryFunc ) ( char* status );
		typedef std::map<std::string, queryFunc> QueryMap;
		typedef bool ( *setFunc ) ( const char* status );
		typedef std::map<std::string, setFunc> SetupMap;
	public:
		GLRenderSetup();

		bool setRenderOption ( const char* name, const char* status );
		bool queryRenderOption ( const char* name, char* status );

	private:
		static bool _setFogEnabled ( const char* status );
		static bool _setFogMode ( const char* status );
		static bool _setFogDensity ( const char* status );
		static bool _setFogStart ( const char* status );
		static bool _setFogEnd ( const char* status );
		static bool _setFogColor ( const char* status );
		static void _queryFogEnabled ( char* status );
		static void _queryFogMode ( char* status );
		static void _queryFogDensity ( char* status );
		static void _queryFogStart ( char* status );
		static void _queryFogEnd ( char* status );
		static void _queryFogColor ( char* status );

	private:
		SetupMap _setupMap;
		QueryMap _queryMap;
	};

	class GLUtils
	{
	public:
		static GLenum getLightEnum ( const int index );
	};
}

#endif
