/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

Implementation by the XKCP contributors, hereby denoted as "the implementer".

For more information, feedback or questions, please refer to the Keccak Team website:
https://keccak.team/

Runtime CPU-feature dispatch for AArch64 (selects the ARMv8.4-A SHA3 backend
when available, otherwise the generic 64-bit backend).

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef _aarch64_dispatch_h_
#define _aarch64_dispatch_h_

void XKCP_EnableAllCpuFeatures(void);
int XKCP_DisableSHA3(void);
int XKCP_ProcessCpuFeatureCommandLineOption(const char *arg);

#endif
