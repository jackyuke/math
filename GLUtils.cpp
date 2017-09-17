#include <unidefs.h>
#include "GLUtils.h"

using namespace RenderEngine;

#define BEGIN_LIGHT_MAP(i) switch((i)) {
#define LIGHT_ENTRY(i) case i: ret = GL_LIGHT##i; break;
#define END_LIGHT_MAP() default: break; }

GLenum GLUtils::getLightEnum ( const int index )
{
	GLenum ret = GL_LIGHT0;

	BEGIN_LIGHT_MAP ( index )
	LIGHT_ENTRY ( 0 )
	LIGHT_ENTRY ( 1 )
	LIGHT_ENTRY ( 2 )
	LIGHT_ENTRY ( 3 )
	LIGHT_ENTRY ( 4 )
	LIGHT_ENTRY ( 5 )
	LIGHT_ENTRY ( 6 )
	LIGHT_ENTRY ( 7 )
	END_LIGHT_MAP()

	return ret;
}

class GLToken
{
public:
	static const char* kON;
	static const char* kOFF;

	//Fog attributes
	static const char* kENABLEFOG;
	static const char* kFOGMODE;
	static const char* kFOGDENSITY;
	static const char* kFOGSTART;
	static const char* kFOGEND;
	static const char* kFOGCOLOR;
	static const char* kEXP;
	static const char* kLINEAR;
	static const char* kEXP2;
};

const char* GLToken::kON = "on";
const char* GLToken::kOFF = "off";
const char* GLToken::kENABLEFOG = "fog enabled";
const char* GLToken::kFOGMODE = "fog mode";
const char* GLToken::kFOGDENSITY = "fog density";
const char* GLToken::kFOGSTART = "fog start";
const char* GLToken::kFOGEND = "fog end";
const char* GLToken::kFOGCOLOR = "fog color";
const char* GLToken::kEXP = "exp";
const char* GLToken::kLINEAR = "linear";
const char* GLToken::kEXP2 = "exp2";

//opengl status query function and map
#define QUERYMAP(name, func) _queryMap[name] = func;

//status setup function and map
#define SETUPMAP(name, func) _setupMap[name] = func;
GLRenderSetup::GLRenderSetup()
{
	// build query map
	QUERYMAP ( GLToken::kENABLEFOG, _queryFogEnabled );
	QUERYMAP ( GLToken::kFOGMODE, _queryFogMode );
	QUERYMAP ( GLToken::kFOGDENSITY, _queryFogDensity );
	QUERYMAP ( GLToken::kFOGSTART, _queryFogStart );
	QUERYMAP ( GLToken::kFOGEND, _queryFogEnd );
	QUERYMAP ( GLToken::kFOGCOLOR, _queryFogColor );
	//build setup map
	SETUPMAP ( GLToken::kENABLEFOG, _setFogEnabled );
	SETUPMAP ( GLToken::kFOGMODE, _setFogMode );
	SETUPMAP ( GLToken::kFOGDENSITY, _setFogDensity );
	SETUPMAP ( GLToken::kFOGSTART, _setFogStart );
	SETUPMAP ( GLToken::kFOGEND, _setFogEnd );
	SETUPMAP ( GLToken::kFOGCOLOR, _setFogColor );
}

bool GLRenderSetup::_setFogEnabled ( const char* status )
{
	if ( strcmp ( status, GLToken::kON ) == 0 )
	{
		glEnable ( GL_FOG );
		return true;
	}
	else if ( strcmp ( status, GLToken::kOFF ) == 0 )
	{
		glDisable ( GL_FOG );
		return true;
	}

	return false;
}

bool GLRenderSetup::_setFogMode ( const char* status )
{
	if ( strcmp ( status, GLToken::kEXP ) == 0 )
		glFogi ( GL_FOG_MODE, GL_EXP );
	else if ( strcmp ( status, GLToken::kEXP2 ) == 0 )
		glFogi ( GL_FOG_MODE, GL_EXP2 );
	else if ( strcmp ( status, GLToken::kLINEAR ) == 0 )
		glFogi ( GL_FOG_MODE, GL_LINEAR );
	else
		return false;

	return true;
}

bool GLRenderSetup::_setFogDensity ( const char* status )
{
	GLfloat density = 0.0f;
	if ( sscanf ( status, "%f", &density ) == 1 )
	{
		if ( density > 1.0f || density < 0.0f )
			return false;

		glFogf ( GL_FOG_DENSITY, density );
		return true;
	}

	return false;
}

bool GLRenderSetup::_setFogStart ( const char* status )
{
	GLfloat start = 0.0f;
	if ( sscanf ( status, "%f", &start ) == 1 )
	{
		glFogf ( GL_FOG_START, start );
		return true;
	}

	return false;
}

bool GLRenderSetup::_setFogEnd ( const char* status )
{
	GLfloat end = 0.0f;
	if ( sscanf ( status, "%f", &end ) == 1 )
	{
		glFogf ( GL_FOG_END, end );
		return true;
	}

	return false;
}

bool GLRenderSetup::_setFogColor ( const char* status )
{
	GLfloat c[4] = {0.0f, 0.0f, 0.0f, 0.0f};

	if ( sscanf ( status, "%f,%f,%f,%f", &c[0], &c[1], &c[2], &c[3] ) == 4 )
	{
		for ( int i = 0; i < 4; i++ )
		{
			if ( c[i] > 1.0f || c[i] < 0.0f )
				return false;
		}
		glFogfv ( GL_FOG_COLOR, c );
		return true;
	}

	return false;
}

bool GLRenderSetup::setRenderOption ( const char* name, const char* status )
{
	if ( name == NULL || status == NULL )
		return false;

	std::string strName = name;
	if ( _setupMap.find ( strName ) != _setupMap.end() )
	{
		setFunc func = _setupMap[strName];
		if ( func )
			return func ( status );
	}
	return false;
}

void GLRenderSetup::_queryFogEnabled ( char* status )
{
	GLboolean fogEnabled = false;
	glGetBooleanv ( GL_FOG, &fogEnabled );
	strcpy ( status, fogEnabled == GL_TRUE ? GLToken::kON : GLToken::kOFF );
}

void GLRenderSetup::_queryFogDensity ( char* status )
{
	GLfloat density = 0.0f;
	glGetFloatv ( GL_FOG_DENSITY, &density );
	sprintf ( status, "%f", density );
}

void GLRenderSetup::_queryFogStart ( char* status )
{
	GLfloat start = 0.0f;
	glGetFloatv ( GL_FOG_START, &start );
	sprintf ( status, "%f", start );
}

void GLRenderSetup::_queryFogEnd ( char* status )
{
	GLfloat end = 0.0f;
	glGetFloatv ( GL_FOG_END, &end );
	sprintf ( status, "%f", end );
}

void GLRenderSetup::_queryFogColor ( char* status )
{
	GLfloat c[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	glGetFloatv ( GL_FOG_COLOR, c );
	sprintf ( status, "%f,%f,%f,%f", c[0], c[1], c[2], c[3] );
}

void GLRenderSetup::_queryFogMode ( char* status )
{
	GLint mode;
	glGetIntegerv ( GL_FOG_MODE, &mode );
	if ( mode == GL_EXP )
		strcpy ( status, GLToken::kEXP );
	else if ( mode == GL_EXP2 )
		strcpy ( status, GLToken::kEXP2 );
	else if ( mode == GL_LINEAR )
		strcpy ( status, GLToken::kLINEAR );
}

bool GLRenderSetup::queryRenderOption ( const char* name, char* status )
{
	if ( name == NULL || status == NULL )
		return false;

	std::string strName = name;
	if ( _queryMap.find ( strName ) != _queryMap.end() )
	{
		queryFunc func = _queryMap[strName];
		if ( func )
		{
			func ( status );
			return true;
		}
	}
	return false;
}
