#include "md2.h"

int md2_init()
{
	// Set up the context
	if (core_init() != RLC_OK)
	{
		core_clean();
		return 1;
	}
	if (pc_param_set_any())
	{
		THROW(ERR_NO_CURVE);
		core_clean();
		return 1;
	}
	return 0;
}
