// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>


#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <windows.h>

#include <vector>
#include <deque>
#include <stack>
#include <queue>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <limits>
#include <algorithm>
#include <functional>
#include <math.h>

#pragma warning(push)
#pragma warning(disable:4201)
#pragma warning(disable:4819)
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#pragma warning(pop)

static const double TOLERANCE = 0.2;
static std::string INPUTPATH;

#include "cgalpackage.h"

//typedef glm::highp_dvec2 Vector2d;
//typedef glm::highp_dvec3 Vector3d;

using namespace std;

// TODO: reference additional headers your program requires here
