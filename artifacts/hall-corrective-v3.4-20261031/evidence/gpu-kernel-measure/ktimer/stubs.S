    .text
    .globl cuArray3DCreate
    .type cuArray3DCreate,@function
cuArray3DCreate:
    jmp *kt_real_cuArray3DCreate(%rip)
    .globl cuArray3DCreate_v2
    .type cuArray3DCreate_v2,@function
cuArray3DCreate_v2:
    jmp *kt_real_cuArray3DCreate_v2(%rip)
    .globl cuArray3DGetDescriptor
    .type cuArray3DGetDescriptor,@function
cuArray3DGetDescriptor:
    jmp *kt_real_cuArray3DGetDescriptor(%rip)
    .globl cuArray3DGetDescriptor_v2
    .type cuArray3DGetDescriptor_v2,@function
cuArray3DGetDescriptor_v2:
    jmp *kt_real_cuArray3DGetDescriptor_v2(%rip)
    .globl cuArrayCreate
    .type cuArrayCreate,@function
cuArrayCreate:
    jmp *kt_real_cuArrayCreate(%rip)
    .globl cuArrayCreate_v2
    .type cuArrayCreate_v2,@function
cuArrayCreate_v2:
    jmp *kt_real_cuArrayCreate_v2(%rip)
    .globl cuArrayDestroy
    .type cuArrayDestroy,@function
cuArrayDestroy:
    jmp *kt_real_cuArrayDestroy(%rip)
    .globl cuArrayGetDescriptor
    .type cuArrayGetDescriptor,@function
cuArrayGetDescriptor:
    jmp *kt_real_cuArrayGetDescriptor(%rip)
    .globl cuArrayGetDescriptor_v2
    .type cuArrayGetDescriptor_v2,@function
cuArrayGetDescriptor_v2:
    jmp *kt_real_cuArrayGetDescriptor_v2(%rip)
    .globl cuArrayGetMemoryRequirements
    .type cuArrayGetMemoryRequirements,@function
cuArrayGetMemoryRequirements:
    jmp *kt_real_cuArrayGetMemoryRequirements(%rip)
    .globl cuArrayGetPlane
    .type cuArrayGetPlane,@function
cuArrayGetPlane:
    jmp *kt_real_cuArrayGetPlane(%rip)
    .globl cuArrayGetSparseProperties
    .type cuArrayGetSparseProperties,@function
cuArrayGetSparseProperties:
    jmp *kt_real_cuArrayGetSparseProperties(%rip)
    .globl cuCheckpointOperationComplete
    .type cuCheckpointOperationComplete,@function
cuCheckpointOperationComplete:
    jmp *kt_real_cuCheckpointOperationComplete(%rip)
    .globl cuCheckpointProcessCheckpoint
    .type cuCheckpointProcessCheckpoint,@function
cuCheckpointProcessCheckpoint:
    jmp *kt_real_cuCheckpointProcessCheckpoint(%rip)
    .globl cuCheckpointProcessGetRestoreThreadId
    .type cuCheckpointProcessGetRestoreThreadId,@function
cuCheckpointProcessGetRestoreThreadId:
    jmp *kt_real_cuCheckpointProcessGetRestoreThreadId(%rip)
    .globl cuCheckpointProcessGetState
    .type cuCheckpointProcessGetState,@function
cuCheckpointProcessGetState:
    jmp *kt_real_cuCheckpointProcessGetState(%rip)
    .globl cuCheckpointProcessLock
    .type cuCheckpointProcessLock,@function
cuCheckpointProcessLock:
    jmp *kt_real_cuCheckpointProcessLock(%rip)
    .globl cuCheckpointProcessRestore
    .type cuCheckpointProcessRestore,@function
cuCheckpointProcessRestore:
    jmp *kt_real_cuCheckpointProcessRestore(%rip)
    .globl cuCheckpointProcessUnlock
    .type cuCheckpointProcessUnlock,@function
cuCheckpointProcessUnlock:
    jmp *kt_real_cuCheckpointProcessUnlock(%rip)
    .globl cuCoredumpDeregisterCompleteCallback
    .type cuCoredumpDeregisterCompleteCallback,@function
cuCoredumpDeregisterCompleteCallback:
    jmp *kt_real_cuCoredumpDeregisterCompleteCallback(%rip)
    .globl cuCoredumpDeregisterStartCallback
    .type cuCoredumpDeregisterStartCallback,@function
cuCoredumpDeregisterStartCallback:
    jmp *kt_real_cuCoredumpDeregisterStartCallback(%rip)
    .globl cuCoredumpGetAttribute
    .type cuCoredumpGetAttribute,@function
cuCoredumpGetAttribute:
    jmp *kt_real_cuCoredumpGetAttribute(%rip)
    .globl cuCoredumpGetAttributeGlobal
    .type cuCoredumpGetAttributeGlobal,@function
cuCoredumpGetAttributeGlobal:
    jmp *kt_real_cuCoredumpGetAttributeGlobal(%rip)
    .globl cuCoredumpRegisterCompleteCallback
    .type cuCoredumpRegisterCompleteCallback,@function
cuCoredumpRegisterCompleteCallback:
    jmp *kt_real_cuCoredumpRegisterCompleteCallback(%rip)
    .globl cuCoredumpRegisterStartCallback
    .type cuCoredumpRegisterStartCallback,@function
cuCoredumpRegisterStartCallback:
    jmp *kt_real_cuCoredumpRegisterStartCallback(%rip)
    .globl cuCoredumpSetAttribute
    .type cuCoredumpSetAttribute,@function
cuCoredumpSetAttribute:
    jmp *kt_real_cuCoredumpSetAttribute(%rip)
    .globl cuCoredumpSetAttributeGlobal
    .type cuCoredumpSetAttributeGlobal,@function
cuCoredumpSetAttributeGlobal:
    jmp *kt_real_cuCoredumpSetAttributeGlobal(%rip)
    .globl cuCtxAttach
    .type cuCtxAttach,@function
cuCtxAttach:
    jmp *kt_real_cuCtxAttach(%rip)
    .globl cuCtxCreate
    .type cuCtxCreate,@function
cuCtxCreate:
    jmp *kt_real_cuCtxCreate(%rip)
    .globl cuCtxCreate_v2
    .type cuCtxCreate_v2,@function
cuCtxCreate_v2:
    jmp *kt_real_cuCtxCreate_v2(%rip)
    .globl cuCtxCreate_v3
    .type cuCtxCreate_v3,@function
cuCtxCreate_v3:
    jmp *kt_real_cuCtxCreate_v3(%rip)
    .globl cuCtxCreate_v4
    .type cuCtxCreate_v4,@function
cuCtxCreate_v4:
    jmp *kt_real_cuCtxCreate_v4(%rip)
    .globl cuCtxDestroy
    .type cuCtxDestroy,@function
cuCtxDestroy:
    jmp *kt_real_cuCtxDestroy(%rip)
    .globl cuCtxDestroy_v2
    .type cuCtxDestroy_v2,@function
cuCtxDestroy_v2:
    jmp *kt_real_cuCtxDestroy_v2(%rip)
    .globl cuCtxDetach
    .type cuCtxDetach,@function
cuCtxDetach:
    jmp *kt_real_cuCtxDetach(%rip)
    .globl cuCtxDisablePeerAccess
    .type cuCtxDisablePeerAccess,@function
cuCtxDisablePeerAccess:
    jmp *kt_real_cuCtxDisablePeerAccess(%rip)
    .globl cuCtxEnablePeerAccess
    .type cuCtxEnablePeerAccess,@function
cuCtxEnablePeerAccess:
    jmp *kt_real_cuCtxEnablePeerAccess(%rip)
    .globl cuCtxFromGreenCtx
    .type cuCtxFromGreenCtx,@function
cuCtxFromGreenCtx:
    jmp *kt_real_cuCtxFromGreenCtx(%rip)
    .globl cuCtxGetApiVersion
    .type cuCtxGetApiVersion,@function
cuCtxGetApiVersion:
    jmp *kt_real_cuCtxGetApiVersion(%rip)
    .globl cuCtxGetCacheConfig
    .type cuCtxGetCacheConfig,@function
cuCtxGetCacheConfig:
    jmp *kt_real_cuCtxGetCacheConfig(%rip)
    .globl cuCtxGetCurrent
    .type cuCtxGetCurrent,@function
cuCtxGetCurrent:
    jmp *kt_real_cuCtxGetCurrent(%rip)
    .globl cuCtxGetDevice
    .type cuCtxGetDevice,@function
cuCtxGetDevice:
    jmp *kt_real_cuCtxGetDevice(%rip)
    .globl cuCtxGetDevice_v2
    .type cuCtxGetDevice_v2,@function
cuCtxGetDevice_v2:
    jmp *kt_real_cuCtxGetDevice_v2(%rip)
    .globl cuCtxGetDevResource
    .type cuCtxGetDevResource,@function
cuCtxGetDevResource:
    jmp *kt_real_cuCtxGetDevResource(%rip)
    .globl cuCtxGetExecAffinity
    .type cuCtxGetExecAffinity,@function
cuCtxGetExecAffinity:
    jmp *kt_real_cuCtxGetExecAffinity(%rip)
    .globl cuCtxGetFlags
    .type cuCtxGetFlags,@function
cuCtxGetFlags:
    jmp *kt_real_cuCtxGetFlags(%rip)
    .globl cuCtxGetId
    .type cuCtxGetId,@function
cuCtxGetId:
    jmp *kt_real_cuCtxGetId(%rip)
    .globl cuCtxGetLimit
    .type cuCtxGetLimit,@function
cuCtxGetLimit:
    jmp *kt_real_cuCtxGetLimit(%rip)
    .globl cuCtxGetSharedMemConfig
    .type cuCtxGetSharedMemConfig,@function
cuCtxGetSharedMemConfig:
    jmp *kt_real_cuCtxGetSharedMemConfig(%rip)
    .globl cuCtxGetStreamPriorityRange
    .type cuCtxGetStreamPriorityRange,@function
cuCtxGetStreamPriorityRange:
    jmp *kt_real_cuCtxGetStreamPriorityRange(%rip)
    .globl cuCtxPopCurrent
    .type cuCtxPopCurrent,@function
cuCtxPopCurrent:
    jmp *kt_real_cuCtxPopCurrent(%rip)
    .globl cuCtxPopCurrent_v2
    .type cuCtxPopCurrent_v2,@function
cuCtxPopCurrent_v2:
    jmp *kt_real_cuCtxPopCurrent_v2(%rip)
    .globl cuCtxPushCurrent
    .type cuCtxPushCurrent,@function
cuCtxPushCurrent:
    jmp *kt_real_cuCtxPushCurrent(%rip)
    .globl cuCtxPushCurrent_v2
    .type cuCtxPushCurrent_v2,@function
cuCtxPushCurrent_v2:
    jmp *kt_real_cuCtxPushCurrent_v2(%rip)
    .globl cuCtxRecordEvent
    .type cuCtxRecordEvent,@function
cuCtxRecordEvent:
    jmp *kt_real_cuCtxRecordEvent(%rip)
    .globl cuCtxResetPersistingL2Cache
    .type cuCtxResetPersistingL2Cache,@function
cuCtxResetPersistingL2Cache:
    jmp *kt_real_cuCtxResetPersistingL2Cache(%rip)
    .globl cuCtxSetCacheConfig
    .type cuCtxSetCacheConfig,@function
cuCtxSetCacheConfig:
    jmp *kt_real_cuCtxSetCacheConfig(%rip)
    .globl cuCtxSetCurrent
    .type cuCtxSetCurrent,@function
cuCtxSetCurrent:
    jmp *kt_real_cuCtxSetCurrent(%rip)
    .globl cuCtxSetFlags
    .type cuCtxSetFlags,@function
cuCtxSetFlags:
    jmp *kt_real_cuCtxSetFlags(%rip)
    .globl cuCtxSetLimit
    .type cuCtxSetLimit,@function
cuCtxSetLimit:
    jmp *kt_real_cuCtxSetLimit(%rip)
    .globl cuCtxSetSharedMemConfig
    .type cuCtxSetSharedMemConfig,@function
cuCtxSetSharedMemConfig:
    jmp *kt_real_cuCtxSetSharedMemConfig(%rip)
    .globl cuCtxSynchronize_v2
    .type cuCtxSynchronize_v2,@function
cuCtxSynchronize_v2:
    jmp *kt_real_cuCtxSynchronize_v2(%rip)
    .globl cuCtxWaitEvent
    .type cuCtxWaitEvent,@function
cuCtxWaitEvent:
    jmp *kt_real_cuCtxWaitEvent(%rip)
    .globl cudbgApiAttach
    .type cudbgApiAttach,@function
cudbgApiAttach:
    jmp *kt_real_cudbgApiAttach(%rip)
    .globl cudbgApiDetach
    .type cudbgApiDetach,@function
cudbgApiDetach:
    jmp *kt_real_cudbgApiDetach(%rip)
    .globl cudbgApiInit
    .type cudbgApiInit,@function
cudbgApiInit:
    jmp *kt_real_cudbgApiInit(%rip)
    .globl cudbgGetAPI
    .type cudbgGetAPI,@function
cudbgGetAPI:
    jmp *kt_real_cudbgGetAPI(%rip)
    .globl cudbgGetAPIVersion
    .type cudbgGetAPIVersion,@function
cudbgGetAPIVersion:
    jmp *kt_real_cudbgGetAPIVersion(%rip)
    .globl cudbgMain
    .type cudbgMain,@function
cudbgMain:
    jmp *kt_real_cudbgMain(%rip)
    .globl cudbgPreInit
    .type cudbgPreInit,@function
cudbgPreInit:
    jmp *kt_real_cudbgPreInit(%rip)
    .globl cudbgReportAttachProcedureFinished
    .type cudbgReportAttachProcedureFinished,@function
cudbgReportAttachProcedureFinished:
    jmp *kt_real_cudbgReportAttachProcedureFinished(%rip)
    .globl cudbgReportDriverApiError
    .type cudbgReportDriverApiError,@function
cudbgReportDriverApiError:
    jmp *kt_real_cudbgReportDriverApiError(%rip)
    .globl cudbgReportDriverInternalError
    .type cudbgReportDriverInternalError,@function
cudbgReportDriverInternalError:
    jmp *kt_real_cudbgReportDriverInternalError(%rip)
    .globl cuDestroyExternalMemory
    .type cuDestroyExternalMemory,@function
cuDestroyExternalMemory:
    jmp *kt_real_cuDestroyExternalMemory(%rip)
    .globl cuDestroyExternalSemaphore
    .type cuDestroyExternalSemaphore,@function
cuDestroyExternalSemaphore:
    jmp *kt_real_cuDestroyExternalSemaphore(%rip)
    .globl cuDeviceCanAccessPeer
    .type cuDeviceCanAccessPeer,@function
cuDeviceCanAccessPeer:
    jmp *kt_real_cuDeviceCanAccessPeer(%rip)
    .globl cuDeviceComputeCapability
    .type cuDeviceComputeCapability,@function
cuDeviceComputeCapability:
    jmp *kt_real_cuDeviceComputeCapability(%rip)
    .globl cuDeviceGet
    .type cuDeviceGet,@function
cuDeviceGet:
    jmp *kt_real_cuDeviceGet(%rip)
    .globl cuDeviceGetAttribute
    .type cuDeviceGetAttribute,@function
cuDeviceGetAttribute:
    jmp *kt_real_cuDeviceGetAttribute(%rip)
    .globl cuDeviceGetByPCIBusId
    .type cuDeviceGetByPCIBusId,@function
cuDeviceGetByPCIBusId:
    jmp *kt_real_cuDeviceGetByPCIBusId(%rip)
    .globl cuDeviceGetCliqueCount
    .type cuDeviceGetCliqueCount,@function
cuDeviceGetCliqueCount:
    jmp *kt_real_cuDeviceGetCliqueCount(%rip)
    .globl cuDeviceGetCliqueInfo
    .type cuDeviceGetCliqueInfo,@function
cuDeviceGetCliqueInfo:
    jmp *kt_real_cuDeviceGetCliqueInfo(%rip)
    .globl cuDeviceGetCount
    .type cuDeviceGetCount,@function
cuDeviceGetCount:
    jmp *kt_real_cuDeviceGetCount(%rip)
    .globl cuDeviceGetDefaultMemPool
    .type cuDeviceGetDefaultMemPool,@function
cuDeviceGetDefaultMemPool:
    jmp *kt_real_cuDeviceGetDefaultMemPool(%rip)
    .globl cuDeviceGetDevResource
    .type cuDeviceGetDevResource,@function
cuDeviceGetDevResource:
    jmp *kt_real_cuDeviceGetDevResource(%rip)
    .globl cuDeviceGetExecAffinitySupport
    .type cuDeviceGetExecAffinitySupport,@function
cuDeviceGetExecAffinitySupport:
    jmp *kt_real_cuDeviceGetExecAffinitySupport(%rip)
    .globl cuDeviceGetFabricClusterUuid
    .type cuDeviceGetFabricClusterUuid,@function
cuDeviceGetFabricClusterUuid:
    jmp *kt_real_cuDeviceGetFabricClusterUuid(%rip)
    .globl cuDeviceGetGraphMemAttribute
    .type cuDeviceGetGraphMemAttribute,@function
cuDeviceGetGraphMemAttribute:
    jmp *kt_real_cuDeviceGetGraphMemAttribute(%rip)
    .globl cuDeviceGetHostAtomicCapabilities
    .type cuDeviceGetHostAtomicCapabilities,@function
cuDeviceGetHostAtomicCapabilities:
    jmp *kt_real_cuDeviceGetHostAtomicCapabilities(%rip)
    .globl cuDeviceGetLuid
    .type cuDeviceGetLuid,@function
cuDeviceGetLuid:
    jmp *kt_real_cuDeviceGetLuid(%rip)
    .globl cuDeviceGetMemPool
    .type cuDeviceGetMemPool,@function
cuDeviceGetMemPool:
    jmp *kt_real_cuDeviceGetMemPool(%rip)
    .globl cuDeviceGetName
    .type cuDeviceGetName,@function
cuDeviceGetName:
    jmp *kt_real_cuDeviceGetName(%rip)
    .globl cuDeviceGetNvSciSyncAttributes
    .type cuDeviceGetNvSciSyncAttributes,@function
cuDeviceGetNvSciSyncAttributes:
    jmp *kt_real_cuDeviceGetNvSciSyncAttributes(%rip)
    .globl cuDeviceGetP2PAtomicCapabilities
    .type cuDeviceGetP2PAtomicCapabilities,@function
cuDeviceGetP2PAtomicCapabilities:
    jmp *kt_real_cuDeviceGetP2PAtomicCapabilities(%rip)
    .globl cuDeviceGetP2PAttribute
    .type cuDeviceGetP2PAttribute,@function
cuDeviceGetP2PAttribute:
    jmp *kt_real_cuDeviceGetP2PAttribute(%rip)
    .globl cuDeviceGetPCIBusId
    .type cuDeviceGetPCIBusId,@function
cuDeviceGetPCIBusId:
    jmp *kt_real_cuDeviceGetPCIBusId(%rip)
    .globl cuDeviceGetProperties
    .type cuDeviceGetProperties,@function
cuDeviceGetProperties:
    jmp *kt_real_cuDeviceGetProperties(%rip)
    .globl cuDeviceGetTexture1DLinearMaxWidth
    .type cuDeviceGetTexture1DLinearMaxWidth,@function
cuDeviceGetTexture1DLinearMaxWidth:
    jmp *kt_real_cuDeviceGetTexture1DLinearMaxWidth(%rip)
    .globl cuDeviceGetUuid
    .type cuDeviceGetUuid,@function
cuDeviceGetUuid:
    jmp *kt_real_cuDeviceGetUuid(%rip)
    .globl cuDeviceGetUuid_v2
    .type cuDeviceGetUuid_v2,@function
cuDeviceGetUuid_v2:
    jmp *kt_real_cuDeviceGetUuid_v2(%rip)
    .globl cuDeviceGraphMemTrim
    .type cuDeviceGraphMemTrim,@function
cuDeviceGraphMemTrim:
    jmp *kt_real_cuDeviceGraphMemTrim(%rip)
    .globl cuDevicePrimaryCtxGetState
    .type cuDevicePrimaryCtxGetState,@function
cuDevicePrimaryCtxGetState:
    jmp *kt_real_cuDevicePrimaryCtxGetState(%rip)
    .globl cuDevicePrimaryCtxRelease
    .type cuDevicePrimaryCtxRelease,@function
cuDevicePrimaryCtxRelease:
    jmp *kt_real_cuDevicePrimaryCtxRelease(%rip)
    .globl cuDevicePrimaryCtxRelease_v2
    .type cuDevicePrimaryCtxRelease_v2,@function
cuDevicePrimaryCtxRelease_v2:
    jmp *kt_real_cuDevicePrimaryCtxRelease_v2(%rip)
    .globl cuDevicePrimaryCtxReset
    .type cuDevicePrimaryCtxReset,@function
cuDevicePrimaryCtxReset:
    jmp *kt_real_cuDevicePrimaryCtxReset(%rip)
    .globl cuDevicePrimaryCtxReset_v2
    .type cuDevicePrimaryCtxReset_v2,@function
cuDevicePrimaryCtxReset_v2:
    jmp *kt_real_cuDevicePrimaryCtxReset_v2(%rip)
    .globl cuDevicePrimaryCtxRetain
    .type cuDevicePrimaryCtxRetain,@function
cuDevicePrimaryCtxRetain:
    jmp *kt_real_cuDevicePrimaryCtxRetain(%rip)
    .globl cuDevicePrimaryCtxSetFlags
    .type cuDevicePrimaryCtxSetFlags,@function
cuDevicePrimaryCtxSetFlags:
    jmp *kt_real_cuDevicePrimaryCtxSetFlags(%rip)
    .globl cuDevicePrimaryCtxSetFlags_v2
    .type cuDevicePrimaryCtxSetFlags_v2,@function
cuDevicePrimaryCtxSetFlags_v2:
    jmp *kt_real_cuDevicePrimaryCtxSetFlags_v2(%rip)
    .globl cuDeviceRegisterAsyncNotification
    .type cuDeviceRegisterAsyncNotification,@function
cuDeviceRegisterAsyncNotification:
    jmp *kt_real_cuDeviceRegisterAsyncNotification(%rip)
    .globl cuDeviceSetGraphMemAttribute
    .type cuDeviceSetGraphMemAttribute,@function
cuDeviceSetGraphMemAttribute:
    jmp *kt_real_cuDeviceSetGraphMemAttribute(%rip)
    .globl cuDeviceSetMemPool
    .type cuDeviceSetMemPool,@function
cuDeviceSetMemPool:
    jmp *kt_real_cuDeviceSetMemPool(%rip)
    .globl cuDeviceTotalMem
    .type cuDeviceTotalMem,@function
cuDeviceTotalMem:
    jmp *kt_real_cuDeviceTotalMem(%rip)
    .globl cuDeviceTotalMem_v2
    .type cuDeviceTotalMem_v2,@function
cuDeviceTotalMem_v2:
    jmp *kt_real_cuDeviceTotalMem_v2(%rip)
    .globl cuDeviceUnregisterAsyncNotification
    .type cuDeviceUnregisterAsyncNotification,@function
cuDeviceUnregisterAsyncNotification:
    jmp *kt_real_cuDeviceUnregisterAsyncNotification(%rip)
    .globl cuDevResourceGenerateDesc
    .type cuDevResourceGenerateDesc,@function
cuDevResourceGenerateDesc:
    jmp *kt_real_cuDevResourceGenerateDesc(%rip)
    .globl cuDevSmResourceSplit
    .type cuDevSmResourceSplit,@function
cuDevSmResourceSplit:
    jmp *kt_real_cuDevSmResourceSplit(%rip)
    .globl cuDevSmResourceSplitByCount
    .type cuDevSmResourceSplitByCount,@function
cuDevSmResourceSplitByCount:
    jmp *kt_real_cuDevSmResourceSplitByCount(%rip)
    .globl cuDriverGetGpuCodeIsaVersion
    .type cuDriverGetGpuCodeIsaVersion,@function
cuDriverGetGpuCodeIsaVersion:
    jmp *kt_real_cuDriverGetGpuCodeIsaVersion(%rip)
    .globl cuDriverGetVersion
    .type cuDriverGetVersion,@function
cuDriverGetVersion:
    jmp *kt_real_cuDriverGetVersion(%rip)
    .globl cuEGLApiInit
    .type cuEGLApiInit,@function
cuEGLApiInit:
    jmp *kt_real_cuEGLApiInit(%rip)
    .globl cuEGLStreamConsumerAcquireFrame
    .type cuEGLStreamConsumerAcquireFrame,@function
cuEGLStreamConsumerAcquireFrame:
    jmp *kt_real_cuEGLStreamConsumerAcquireFrame(%rip)
    .globl cuEGLStreamConsumerConnect
    .type cuEGLStreamConsumerConnect,@function
cuEGLStreamConsumerConnect:
    jmp *kt_real_cuEGLStreamConsumerConnect(%rip)
    .globl cuEGLStreamConsumerConnectWithFlags
    .type cuEGLStreamConsumerConnectWithFlags,@function
cuEGLStreamConsumerConnectWithFlags:
    jmp *kt_real_cuEGLStreamConsumerConnectWithFlags(%rip)
    .globl cuEGLStreamConsumerDisconnect
    .type cuEGLStreamConsumerDisconnect,@function
cuEGLStreamConsumerDisconnect:
    jmp *kt_real_cuEGLStreamConsumerDisconnect(%rip)
    .globl cuEGLStreamConsumerReleaseFrame
    .type cuEGLStreamConsumerReleaseFrame,@function
cuEGLStreamConsumerReleaseFrame:
    jmp *kt_real_cuEGLStreamConsumerReleaseFrame(%rip)
    .globl cuEGLStreamProducerConnect
    .type cuEGLStreamProducerConnect,@function
cuEGLStreamProducerConnect:
    jmp *kt_real_cuEGLStreamProducerConnect(%rip)
    .globl cuEGLStreamProducerDisconnect
    .type cuEGLStreamProducerDisconnect,@function
cuEGLStreamProducerDisconnect:
    jmp *kt_real_cuEGLStreamProducerDisconnect(%rip)
    .globl cuEGLStreamProducerPresentFrame
    .type cuEGLStreamProducerPresentFrame,@function
cuEGLStreamProducerPresentFrame:
    jmp *kt_real_cuEGLStreamProducerPresentFrame(%rip)
    .globl cuEGLStreamProducerReturnFrame
    .type cuEGLStreamProducerReturnFrame,@function
cuEGLStreamProducerReturnFrame:
    jmp *kt_real_cuEGLStreamProducerReturnFrame(%rip)
    .globl cuEventCreate
    .type cuEventCreate,@function
cuEventCreate:
    jmp *kt_real_cuEventCreate(%rip)
    .globl cuEventDestroy
    .type cuEventDestroy,@function
cuEventDestroy:
    jmp *kt_real_cuEventDestroy(%rip)
    .globl cuEventDestroy_v2
    .type cuEventDestroy_v2,@function
cuEventDestroy_v2:
    jmp *kt_real_cuEventDestroy_v2(%rip)
    .globl cuEventElapsedTime
    .type cuEventElapsedTime,@function
cuEventElapsedTime:
    jmp *kt_real_cuEventElapsedTime(%rip)
    .globl cuEventElapsedTime_v2
    .type cuEventElapsedTime_v2,@function
cuEventElapsedTime_v2:
    jmp *kt_real_cuEventElapsedTime_v2(%rip)
    .globl cuEventQuery
    .type cuEventQuery,@function
cuEventQuery:
    jmp *kt_real_cuEventQuery(%rip)
    .globl cuEventRecord
    .type cuEventRecord,@function
cuEventRecord:
    jmp *kt_real_cuEventRecord(%rip)
    .globl cuEventRecord_ptsz
    .type cuEventRecord_ptsz,@function
cuEventRecord_ptsz:
    jmp *kt_real_cuEventRecord_ptsz(%rip)
    .globl cuEventRecordWithFlags
    .type cuEventRecordWithFlags,@function
cuEventRecordWithFlags:
    jmp *kt_real_cuEventRecordWithFlags(%rip)
    .globl cuEventRecordWithFlags_ptsz
    .type cuEventRecordWithFlags_ptsz,@function
cuEventRecordWithFlags_ptsz:
    jmp *kt_real_cuEventRecordWithFlags_ptsz(%rip)
    .globl cuEventSynchronize
    .type cuEventSynchronize,@function
cuEventSynchronize:
    jmp *kt_real_cuEventSynchronize(%rip)
    .globl cuExternalMemoryGetMappedBuffer
    .type cuExternalMemoryGetMappedBuffer,@function
cuExternalMemoryGetMappedBuffer:
    jmp *kt_real_cuExternalMemoryGetMappedBuffer(%rip)
    .globl cuExternalMemoryGetMappedMipmappedArray
    .type cuExternalMemoryGetMappedMipmappedArray,@function
cuExternalMemoryGetMappedMipmappedArray:
    jmp *kt_real_cuExternalMemoryGetMappedMipmappedArray(%rip)
    .globl cuFlushGPUDirectRDMAWrites
    .type cuFlushGPUDirectRDMAWrites,@function
cuFlushGPUDirectRDMAWrites:
    jmp *kt_real_cuFlushGPUDirectRDMAWrites(%rip)
    .globl cuFuncGetAttribute
    .type cuFuncGetAttribute,@function
cuFuncGetAttribute:
    jmp *kt_real_cuFuncGetAttribute(%rip)
    .globl cuFuncGetModule
    .type cuFuncGetModule,@function
cuFuncGetModule:
    jmp *kt_real_cuFuncGetModule(%rip)
    .globl cuFuncGetName
    .type cuFuncGetName,@function
cuFuncGetName:
    jmp *kt_real_cuFuncGetName(%rip)
    .globl cuFuncGetParamCount
    .type cuFuncGetParamCount,@function
cuFuncGetParamCount:
    jmp *kt_real_cuFuncGetParamCount(%rip)
    .globl cuFuncGetParamInfo
    .type cuFuncGetParamInfo,@function
cuFuncGetParamInfo:
    jmp *kt_real_cuFuncGetParamInfo(%rip)
    .globl cuFuncIsLoaded
    .type cuFuncIsLoaded,@function
cuFuncIsLoaded:
    jmp *kt_real_cuFuncIsLoaded(%rip)
    .globl cuFuncLoad
    .type cuFuncLoad,@function
cuFuncLoad:
    jmp *kt_real_cuFuncLoad(%rip)
    .globl cuFuncSetAttribute
    .type cuFuncSetAttribute,@function
cuFuncSetAttribute:
    jmp *kt_real_cuFuncSetAttribute(%rip)
    .globl cuFuncSetBlockShape
    .type cuFuncSetBlockShape,@function
cuFuncSetBlockShape:
    jmp *kt_real_cuFuncSetBlockShape(%rip)
    .globl cuFuncSetCacheConfig
    .type cuFuncSetCacheConfig,@function
cuFuncSetCacheConfig:
    jmp *kt_real_cuFuncSetCacheConfig(%rip)
    .globl cuFuncSetSharedMemConfig
    .type cuFuncSetSharedMemConfig,@function
cuFuncSetSharedMemConfig:
    jmp *kt_real_cuFuncSetSharedMemConfig(%rip)
    .globl cuFuncSetSharedSize
    .type cuFuncSetSharedSize,@function
cuFuncSetSharedSize:
    jmp *kt_real_cuFuncSetSharedSize(%rip)
    .globl cuGetErrorName
    .type cuGetErrorName,@function
cuGetErrorName:
    jmp *kt_real_cuGetErrorName(%rip)
    .globl cuGetErrorString
    .type cuGetErrorString,@function
cuGetErrorString:
    jmp *kt_real_cuGetErrorString(%rip)
    .globl cuGetExportTable
    .type cuGetExportTable,@function
cuGetExportTable:
    jmp *kt_real_cuGetExportTable(%rip)
    .globl cuGetProcAddress
    .type cuGetProcAddress,@function
cuGetProcAddress:
    jmp *kt_real_cuGetProcAddress(%rip)
    .globl cuGetProcAddress_v2
    .type cuGetProcAddress_v2,@function
cuGetProcAddress_v2:
    jmp *kt_real_cuGetProcAddress_v2(%rip)
    .globl cuGLCtxCreate
    .type cuGLCtxCreate,@function
cuGLCtxCreate:
    jmp *kt_real_cuGLCtxCreate(%rip)
    .globl cuGLCtxCreate_v2
    .type cuGLCtxCreate_v2,@function
cuGLCtxCreate_v2:
    jmp *kt_real_cuGLCtxCreate_v2(%rip)
    .globl cuGLGetDevices
    .type cuGLGetDevices,@function
cuGLGetDevices:
    jmp *kt_real_cuGLGetDevices(%rip)
    .globl cuGLGetDevices_v2
    .type cuGLGetDevices_v2,@function
cuGLGetDevices_v2:
    jmp *kt_real_cuGLGetDevices_v2(%rip)
    .globl cuGLInit
    .type cuGLInit,@function
cuGLInit:
    jmp *kt_real_cuGLInit(%rip)
    .globl cuGLMapBufferObject
    .type cuGLMapBufferObject,@function
cuGLMapBufferObject:
    jmp *kt_real_cuGLMapBufferObject(%rip)
    .globl cuGLMapBufferObjectAsync
    .type cuGLMapBufferObjectAsync,@function
cuGLMapBufferObjectAsync:
    jmp *kt_real_cuGLMapBufferObjectAsync(%rip)
    .globl cuGLMapBufferObjectAsync_v2
    .type cuGLMapBufferObjectAsync_v2,@function
cuGLMapBufferObjectAsync_v2:
    jmp *kt_real_cuGLMapBufferObjectAsync_v2(%rip)
    .globl cuGLMapBufferObjectAsync_v2_ptsz
    .type cuGLMapBufferObjectAsync_v2_ptsz,@function
cuGLMapBufferObjectAsync_v2_ptsz:
    jmp *kt_real_cuGLMapBufferObjectAsync_v2_ptsz(%rip)
    .globl cuGLMapBufferObject_v2
    .type cuGLMapBufferObject_v2,@function
cuGLMapBufferObject_v2:
    jmp *kt_real_cuGLMapBufferObject_v2(%rip)
    .globl cuGLMapBufferObject_v2_ptds
    .type cuGLMapBufferObject_v2_ptds,@function
cuGLMapBufferObject_v2_ptds:
    jmp *kt_real_cuGLMapBufferObject_v2_ptds(%rip)
    .globl cuGLRegisterBufferObject
    .type cuGLRegisterBufferObject,@function
cuGLRegisterBufferObject:
    jmp *kt_real_cuGLRegisterBufferObject(%rip)
    .globl cuGLSetBufferObjectMapFlags
    .type cuGLSetBufferObjectMapFlags,@function
cuGLSetBufferObjectMapFlags:
    jmp *kt_real_cuGLSetBufferObjectMapFlags(%rip)
    .globl cuGLUnmapBufferObject
    .type cuGLUnmapBufferObject,@function
cuGLUnmapBufferObject:
    jmp *kt_real_cuGLUnmapBufferObject(%rip)
    .globl cuGLUnmapBufferObjectAsync
    .type cuGLUnmapBufferObjectAsync,@function
cuGLUnmapBufferObjectAsync:
    jmp *kt_real_cuGLUnmapBufferObjectAsync(%rip)
    .globl cuGLUnregisterBufferObject
    .type cuGLUnregisterBufferObject,@function
cuGLUnregisterBufferObject:
    jmp *kt_real_cuGLUnregisterBufferObject(%rip)
    .globl cuGraphAddBatchMemOpNode
    .type cuGraphAddBatchMemOpNode,@function
cuGraphAddBatchMemOpNode:
    jmp *kt_real_cuGraphAddBatchMemOpNode(%rip)
    .globl cuGraphAddChildGraphNode
    .type cuGraphAddChildGraphNode,@function
cuGraphAddChildGraphNode:
    jmp *kt_real_cuGraphAddChildGraphNode(%rip)
    .globl cuGraphAddDependencies
    .type cuGraphAddDependencies,@function
cuGraphAddDependencies:
    jmp *kt_real_cuGraphAddDependencies(%rip)
    .globl cuGraphAddDependencies_v2
    .type cuGraphAddDependencies_v2,@function
cuGraphAddDependencies_v2:
    jmp *kt_real_cuGraphAddDependencies_v2(%rip)
    .globl cuGraphAddEmptyNode
    .type cuGraphAddEmptyNode,@function
cuGraphAddEmptyNode:
    jmp *kt_real_cuGraphAddEmptyNode(%rip)
    .globl cuGraphAddEventRecordNode
    .type cuGraphAddEventRecordNode,@function
cuGraphAddEventRecordNode:
    jmp *kt_real_cuGraphAddEventRecordNode(%rip)
    .globl cuGraphAddEventWaitNode
    .type cuGraphAddEventWaitNode,@function
cuGraphAddEventWaitNode:
    jmp *kt_real_cuGraphAddEventWaitNode(%rip)
    .globl cuGraphAddExternalSemaphoresSignalNode
    .type cuGraphAddExternalSemaphoresSignalNode,@function
cuGraphAddExternalSemaphoresSignalNode:
    jmp *kt_real_cuGraphAddExternalSemaphoresSignalNode(%rip)
    .globl cuGraphAddExternalSemaphoresWaitNode
    .type cuGraphAddExternalSemaphoresWaitNode,@function
cuGraphAddExternalSemaphoresWaitNode:
    jmp *kt_real_cuGraphAddExternalSemaphoresWaitNode(%rip)
    .globl cuGraphAddHostNode
    .type cuGraphAddHostNode,@function
cuGraphAddHostNode:
    jmp *kt_real_cuGraphAddHostNode(%rip)
    .globl cuGraphAddKernelNode
    .type cuGraphAddKernelNode,@function
cuGraphAddKernelNode:
    jmp *kt_real_cuGraphAddKernelNode(%rip)
    .globl cuGraphAddKernelNode_v2
    .type cuGraphAddKernelNode_v2,@function
cuGraphAddKernelNode_v2:
    jmp *kt_real_cuGraphAddKernelNode_v2(%rip)
    .globl cuGraphAddMemAllocNode
    .type cuGraphAddMemAllocNode,@function
cuGraphAddMemAllocNode:
    jmp *kt_real_cuGraphAddMemAllocNode(%rip)
    .globl cuGraphAddMemcpyNode
    .type cuGraphAddMemcpyNode,@function
cuGraphAddMemcpyNode:
    jmp *kt_real_cuGraphAddMemcpyNode(%rip)
    .globl cuGraphAddMemFreeNode
    .type cuGraphAddMemFreeNode,@function
cuGraphAddMemFreeNode:
    jmp *kt_real_cuGraphAddMemFreeNode(%rip)
    .globl cuGraphAddMemsetNode
    .type cuGraphAddMemsetNode,@function
cuGraphAddMemsetNode:
    jmp *kt_real_cuGraphAddMemsetNode(%rip)
    .globl cuGraphAddNode
    .type cuGraphAddNode,@function
cuGraphAddNode:
    jmp *kt_real_cuGraphAddNode(%rip)
    .globl cuGraphAddNode_v2
    .type cuGraphAddNode_v2,@function
cuGraphAddNode_v2:
    jmp *kt_real_cuGraphAddNode_v2(%rip)
    .globl cuGraphAddNode_v3
    .type cuGraphAddNode_v3,@function
cuGraphAddNode_v3:
    jmp *kt_real_cuGraphAddNode_v3(%rip)
    .globl cuGraphBatchMemOpNodeGetParams
    .type cuGraphBatchMemOpNodeGetParams,@function
cuGraphBatchMemOpNodeGetParams:
    jmp *kt_real_cuGraphBatchMemOpNodeGetParams(%rip)
    .globl cuGraphBatchMemOpNodeSetParams
    .type cuGraphBatchMemOpNodeSetParams,@function
cuGraphBatchMemOpNodeSetParams:
    jmp *kt_real_cuGraphBatchMemOpNodeSetParams(%rip)
    .globl cuGraphChildGraphNodeGetGraph
    .type cuGraphChildGraphNodeGetGraph,@function
cuGraphChildGraphNodeGetGraph:
    jmp *kt_real_cuGraphChildGraphNodeGetGraph(%rip)
    .globl cuGraphClone
    .type cuGraphClone,@function
cuGraphClone:
    jmp *kt_real_cuGraphClone(%rip)
    .globl cuGraphConditionalHandleCreate
    .type cuGraphConditionalHandleCreate,@function
cuGraphConditionalHandleCreate:
    jmp *kt_real_cuGraphConditionalHandleCreate(%rip)
    .globl cuGraphCreate
    .type cuGraphCreate,@function
cuGraphCreate:
    jmp *kt_real_cuGraphCreate(%rip)
    .globl cuGraphDebugDotPrint
    .type cuGraphDebugDotPrint,@function
cuGraphDebugDotPrint:
    jmp *kt_real_cuGraphDebugDotPrint(%rip)
    .globl cuGraphDestroy
    .type cuGraphDestroy,@function
cuGraphDestroy:
    jmp *kt_real_cuGraphDestroy(%rip)
    .globl cuGraphDestroyNode
    .type cuGraphDestroyNode,@function
cuGraphDestroyNode:
    jmp *kt_real_cuGraphDestroyNode(%rip)
    .globl cuGraphEventRecordNodeGetEvent
    .type cuGraphEventRecordNodeGetEvent,@function
cuGraphEventRecordNodeGetEvent:
    jmp *kt_real_cuGraphEventRecordNodeGetEvent(%rip)
    .globl cuGraphEventRecordNodeSetEvent
    .type cuGraphEventRecordNodeSetEvent,@function
cuGraphEventRecordNodeSetEvent:
    jmp *kt_real_cuGraphEventRecordNodeSetEvent(%rip)
    .globl cuGraphEventWaitNodeGetEvent
    .type cuGraphEventWaitNodeGetEvent,@function
cuGraphEventWaitNodeGetEvent:
    jmp *kt_real_cuGraphEventWaitNodeGetEvent(%rip)
    .globl cuGraphEventWaitNodeSetEvent
    .type cuGraphEventWaitNodeSetEvent,@function
cuGraphEventWaitNodeSetEvent:
    jmp *kt_real_cuGraphEventWaitNodeSetEvent(%rip)
    .globl cuGraphExecBatchMemOpNodeSetParams
    .type cuGraphExecBatchMemOpNodeSetParams,@function
cuGraphExecBatchMemOpNodeSetParams:
    jmp *kt_real_cuGraphExecBatchMemOpNodeSetParams(%rip)
    .globl cuGraphExecChildGraphNodeSetParams
    .type cuGraphExecChildGraphNodeSetParams,@function
cuGraphExecChildGraphNodeSetParams:
    jmp *kt_real_cuGraphExecChildGraphNodeSetParams(%rip)
    .globl cuGraphExecDestroy
    .type cuGraphExecDestroy,@function
cuGraphExecDestroy:
    jmp *kt_real_cuGraphExecDestroy(%rip)
    .globl cuGraphExecEventRecordNodeSetEvent
    .type cuGraphExecEventRecordNodeSetEvent,@function
cuGraphExecEventRecordNodeSetEvent:
    jmp *kt_real_cuGraphExecEventRecordNodeSetEvent(%rip)
    .globl cuGraphExecEventWaitNodeSetEvent
    .type cuGraphExecEventWaitNodeSetEvent,@function
cuGraphExecEventWaitNodeSetEvent:
    jmp *kt_real_cuGraphExecEventWaitNodeSetEvent(%rip)
    .globl cuGraphExecExternalSemaphoresSignalNodeSetParams
    .type cuGraphExecExternalSemaphoresSignalNodeSetParams,@function
cuGraphExecExternalSemaphoresSignalNodeSetParams:
    jmp *kt_real_cuGraphExecExternalSemaphoresSignalNodeSetParams(%rip)
    .globl cuGraphExecExternalSemaphoresWaitNodeSetParams
    .type cuGraphExecExternalSemaphoresWaitNodeSetParams,@function
cuGraphExecExternalSemaphoresWaitNodeSetParams:
    jmp *kt_real_cuGraphExecExternalSemaphoresWaitNodeSetParams(%rip)
    .globl cuGraphExecGetFlags
    .type cuGraphExecGetFlags,@function
cuGraphExecGetFlags:
    jmp *kt_real_cuGraphExecGetFlags(%rip)
    .globl cuGraphExecGetId
    .type cuGraphExecGetId,@function
cuGraphExecGetId:
    jmp *kt_real_cuGraphExecGetId(%rip)
    .globl cuGraphExecHostNodeSetParams
    .type cuGraphExecHostNodeSetParams,@function
cuGraphExecHostNodeSetParams:
    jmp *kt_real_cuGraphExecHostNodeSetParams(%rip)
    .globl cuGraphExecKernelNodeSetParams
    .type cuGraphExecKernelNodeSetParams,@function
cuGraphExecKernelNodeSetParams:
    jmp *kt_real_cuGraphExecKernelNodeSetParams(%rip)
    .globl cuGraphExecKernelNodeSetParams_v2
    .type cuGraphExecKernelNodeSetParams_v2,@function
cuGraphExecKernelNodeSetParams_v2:
    jmp *kt_real_cuGraphExecKernelNodeSetParams_v2(%rip)
    .globl cuGraphExecMemcpyNodeSetParams
    .type cuGraphExecMemcpyNodeSetParams,@function
cuGraphExecMemcpyNodeSetParams:
    jmp *kt_real_cuGraphExecMemcpyNodeSetParams(%rip)
    .globl cuGraphExecMemsetNodeSetParams
    .type cuGraphExecMemsetNodeSetParams,@function
cuGraphExecMemsetNodeSetParams:
    jmp *kt_real_cuGraphExecMemsetNodeSetParams(%rip)
    .globl cuGraphExecNodeSetParams
    .type cuGraphExecNodeSetParams,@function
cuGraphExecNodeSetParams:
    jmp *kt_real_cuGraphExecNodeSetParams(%rip)
    .globl cuGraphExecUpdate
    .type cuGraphExecUpdate,@function
cuGraphExecUpdate:
    jmp *kt_real_cuGraphExecUpdate(%rip)
    .globl cuGraphExecUpdate_v2
    .type cuGraphExecUpdate_v2,@function
cuGraphExecUpdate_v2:
    jmp *kt_real_cuGraphExecUpdate_v2(%rip)
    .globl cuGraphExternalSemaphoresSignalNodeGetParams
    .type cuGraphExternalSemaphoresSignalNodeGetParams,@function
cuGraphExternalSemaphoresSignalNodeGetParams:
    jmp *kt_real_cuGraphExternalSemaphoresSignalNodeGetParams(%rip)
    .globl cuGraphExternalSemaphoresSignalNodeSetParams
    .type cuGraphExternalSemaphoresSignalNodeSetParams,@function
cuGraphExternalSemaphoresSignalNodeSetParams:
    jmp *kt_real_cuGraphExternalSemaphoresSignalNodeSetParams(%rip)
    .globl cuGraphExternalSemaphoresWaitNodeGetParams
    .type cuGraphExternalSemaphoresWaitNodeGetParams,@function
cuGraphExternalSemaphoresWaitNodeGetParams:
    jmp *kt_real_cuGraphExternalSemaphoresWaitNodeGetParams(%rip)
    .globl cuGraphExternalSemaphoresWaitNodeSetParams
    .type cuGraphExternalSemaphoresWaitNodeSetParams,@function
cuGraphExternalSemaphoresWaitNodeSetParams:
    jmp *kt_real_cuGraphExternalSemaphoresWaitNodeSetParams(%rip)
    .globl cuGraphGetEdges
    .type cuGraphGetEdges,@function
cuGraphGetEdges:
    jmp *kt_real_cuGraphGetEdges(%rip)
    .globl cuGraphGetEdges_v2
    .type cuGraphGetEdges_v2,@function
cuGraphGetEdges_v2:
    jmp *kt_real_cuGraphGetEdges_v2(%rip)
    .globl cuGraphGetId
    .type cuGraphGetId,@function
cuGraphGetId:
    jmp *kt_real_cuGraphGetId(%rip)
    .globl cuGraphGetNodes
    .type cuGraphGetNodes,@function
cuGraphGetNodes:
    jmp *kt_real_cuGraphGetNodes(%rip)
    .globl cuGraphGetRootNodes
    .type cuGraphGetRootNodes,@function
cuGraphGetRootNodes:
    jmp *kt_real_cuGraphGetRootNodes(%rip)
    .globl cuGraphHostNodeGetParams
    .type cuGraphHostNodeGetParams,@function
cuGraphHostNodeGetParams:
    jmp *kt_real_cuGraphHostNodeGetParams(%rip)
    .globl cuGraphHostNodeSetParams
    .type cuGraphHostNodeSetParams,@function
cuGraphHostNodeSetParams:
    jmp *kt_real_cuGraphHostNodeSetParams(%rip)
    .globl cuGraphicsEGLRegisterImage
    .type cuGraphicsEGLRegisterImage,@function
cuGraphicsEGLRegisterImage:
    jmp *kt_real_cuGraphicsEGLRegisterImage(%rip)
    .globl cuGraphicsGLRegisterBuffer
    .type cuGraphicsGLRegisterBuffer,@function
cuGraphicsGLRegisterBuffer:
    jmp *kt_real_cuGraphicsGLRegisterBuffer(%rip)
    .globl cuGraphicsGLRegisterImage
    .type cuGraphicsGLRegisterImage,@function
cuGraphicsGLRegisterImage:
    jmp *kt_real_cuGraphicsGLRegisterImage(%rip)
    .globl cuGraphicsMapResources
    .type cuGraphicsMapResources,@function
cuGraphicsMapResources:
    jmp *kt_real_cuGraphicsMapResources(%rip)
    .globl cuGraphicsMapResources_ptsz
    .type cuGraphicsMapResources_ptsz,@function
cuGraphicsMapResources_ptsz:
    jmp *kt_real_cuGraphicsMapResources_ptsz(%rip)
    .globl cuGraphicsResourceGetMappedEglFrame
    .type cuGraphicsResourceGetMappedEglFrame,@function
cuGraphicsResourceGetMappedEglFrame:
    jmp *kt_real_cuGraphicsResourceGetMappedEglFrame(%rip)
    .globl cuGraphicsResourceGetMappedMipmappedArray
    .type cuGraphicsResourceGetMappedMipmappedArray,@function
cuGraphicsResourceGetMappedMipmappedArray:
    jmp *kt_real_cuGraphicsResourceGetMappedMipmappedArray(%rip)
    .globl cuGraphicsResourceGetMappedPointer
    .type cuGraphicsResourceGetMappedPointer,@function
cuGraphicsResourceGetMappedPointer:
    jmp *kt_real_cuGraphicsResourceGetMappedPointer(%rip)
    .globl cuGraphicsResourceGetMappedPointer_v2
    .type cuGraphicsResourceGetMappedPointer_v2,@function
cuGraphicsResourceGetMappedPointer_v2:
    jmp *kt_real_cuGraphicsResourceGetMappedPointer_v2(%rip)
    .globl cuGraphicsResourceSetMapFlags
    .type cuGraphicsResourceSetMapFlags,@function
cuGraphicsResourceSetMapFlags:
    jmp *kt_real_cuGraphicsResourceSetMapFlags(%rip)
    .globl cuGraphicsResourceSetMapFlags_v2
    .type cuGraphicsResourceSetMapFlags_v2,@function
cuGraphicsResourceSetMapFlags_v2:
    jmp *kt_real_cuGraphicsResourceSetMapFlags_v2(%rip)
    .globl cuGraphicsSubResourceGetMappedArray
    .type cuGraphicsSubResourceGetMappedArray,@function
cuGraphicsSubResourceGetMappedArray:
    jmp *kt_real_cuGraphicsSubResourceGetMappedArray(%rip)
    .globl cuGraphicsUnmapResources
    .type cuGraphicsUnmapResources,@function
cuGraphicsUnmapResources:
    jmp *kt_real_cuGraphicsUnmapResources(%rip)
    .globl cuGraphicsUnmapResources_ptsz
    .type cuGraphicsUnmapResources_ptsz,@function
cuGraphicsUnmapResources_ptsz:
    jmp *kt_real_cuGraphicsUnmapResources_ptsz(%rip)
    .globl cuGraphicsUnregisterResource
    .type cuGraphicsUnregisterResource,@function
cuGraphicsUnregisterResource:
    jmp *kt_real_cuGraphicsUnregisterResource(%rip)
    .globl cuGraphicsVDPAURegisterOutputSurface
    .type cuGraphicsVDPAURegisterOutputSurface,@function
cuGraphicsVDPAURegisterOutputSurface:
    jmp *kt_real_cuGraphicsVDPAURegisterOutputSurface(%rip)
    .globl cuGraphicsVDPAURegisterVideoSurface
    .type cuGraphicsVDPAURegisterVideoSurface,@function
cuGraphicsVDPAURegisterVideoSurface:
    jmp *kt_real_cuGraphicsVDPAURegisterVideoSurface(%rip)
    .globl cuGraphInstantiate
    .type cuGraphInstantiate,@function
cuGraphInstantiate:
    jmp *kt_real_cuGraphInstantiate(%rip)
    .globl cuGraphInstantiate_v2
    .type cuGraphInstantiate_v2,@function
cuGraphInstantiate_v2:
    jmp *kt_real_cuGraphInstantiate_v2(%rip)
    .globl cuGraphInstantiateWithFlags
    .type cuGraphInstantiateWithFlags,@function
cuGraphInstantiateWithFlags:
    jmp *kt_real_cuGraphInstantiateWithFlags(%rip)
    .globl cuGraphInstantiateWithParams
    .type cuGraphInstantiateWithParams,@function
cuGraphInstantiateWithParams:
    jmp *kt_real_cuGraphInstantiateWithParams(%rip)
    .globl cuGraphInstantiateWithParams_ptsz
    .type cuGraphInstantiateWithParams_ptsz,@function
cuGraphInstantiateWithParams_ptsz:
    jmp *kt_real_cuGraphInstantiateWithParams_ptsz(%rip)
    .globl cuGraphKernelNodeCopyAttributes
    .type cuGraphKernelNodeCopyAttributes,@function
cuGraphKernelNodeCopyAttributes:
    jmp *kt_real_cuGraphKernelNodeCopyAttributes(%rip)
    .globl cuGraphKernelNodeGetAttribute
    .type cuGraphKernelNodeGetAttribute,@function
cuGraphKernelNodeGetAttribute:
    jmp *kt_real_cuGraphKernelNodeGetAttribute(%rip)
    .globl cuGraphKernelNodeGetParams
    .type cuGraphKernelNodeGetParams,@function
cuGraphKernelNodeGetParams:
    jmp *kt_real_cuGraphKernelNodeGetParams(%rip)
    .globl cuGraphKernelNodeGetParams_v2
    .type cuGraphKernelNodeGetParams_v2,@function
cuGraphKernelNodeGetParams_v2:
    jmp *kt_real_cuGraphKernelNodeGetParams_v2(%rip)
    .globl cuGraphKernelNodeSetAttribute
    .type cuGraphKernelNodeSetAttribute,@function
cuGraphKernelNodeSetAttribute:
    jmp *kt_real_cuGraphKernelNodeSetAttribute(%rip)
    .globl cuGraphKernelNodeSetParams
    .type cuGraphKernelNodeSetParams,@function
cuGraphKernelNodeSetParams:
    jmp *kt_real_cuGraphKernelNodeSetParams(%rip)
    .globl cuGraphKernelNodeSetParams_v2
    .type cuGraphKernelNodeSetParams_v2,@function
cuGraphKernelNodeSetParams_v2:
    jmp *kt_real_cuGraphKernelNodeSetParams_v2(%rip)
    .globl cuGraphLaunch
    .type cuGraphLaunch,@function
cuGraphLaunch:
    jmp *kt_real_cuGraphLaunch(%rip)
    .globl cuGraphLaunch_ptsz
    .type cuGraphLaunch_ptsz,@function
cuGraphLaunch_ptsz:
    jmp *kt_real_cuGraphLaunch_ptsz(%rip)
    .globl cuGraphMemAllocNodeGetParams
    .type cuGraphMemAllocNodeGetParams,@function
cuGraphMemAllocNodeGetParams:
    jmp *kt_real_cuGraphMemAllocNodeGetParams(%rip)
    .globl cuGraphMemcpyNodeGetParams
    .type cuGraphMemcpyNodeGetParams,@function
cuGraphMemcpyNodeGetParams:
    jmp *kt_real_cuGraphMemcpyNodeGetParams(%rip)
    .globl cuGraphMemcpyNodeSetParams
    .type cuGraphMemcpyNodeSetParams,@function
cuGraphMemcpyNodeSetParams:
    jmp *kt_real_cuGraphMemcpyNodeSetParams(%rip)
    .globl cuGraphMemFreeNodeGetParams
    .type cuGraphMemFreeNodeGetParams,@function
cuGraphMemFreeNodeGetParams:
    jmp *kt_real_cuGraphMemFreeNodeGetParams(%rip)
    .globl cuGraphMemsetNodeGetParams
    .type cuGraphMemsetNodeGetParams,@function
cuGraphMemsetNodeGetParams:
    jmp *kt_real_cuGraphMemsetNodeGetParams(%rip)
    .globl cuGraphMemsetNodeSetParams
    .type cuGraphMemsetNodeSetParams,@function
cuGraphMemsetNodeSetParams:
    jmp *kt_real_cuGraphMemsetNodeSetParams(%rip)
    .globl cuGraphNodeFindInClone
    .type cuGraphNodeFindInClone,@function
cuGraphNodeFindInClone:
    jmp *kt_real_cuGraphNodeFindInClone(%rip)
    .globl cuGraphNodeGetContainingGraph
    .type cuGraphNodeGetContainingGraph,@function
cuGraphNodeGetContainingGraph:
    jmp *kt_real_cuGraphNodeGetContainingGraph(%rip)
    .globl cuGraphNodeGetDependencies
    .type cuGraphNodeGetDependencies,@function
cuGraphNodeGetDependencies:
    jmp *kt_real_cuGraphNodeGetDependencies(%rip)
    .globl cuGraphNodeGetDependencies_v2
    .type cuGraphNodeGetDependencies_v2,@function
cuGraphNodeGetDependencies_v2:
    jmp *kt_real_cuGraphNodeGetDependencies_v2(%rip)
    .globl cuGraphNodeGetDependentNodes
    .type cuGraphNodeGetDependentNodes,@function
cuGraphNodeGetDependentNodes:
    jmp *kt_real_cuGraphNodeGetDependentNodes(%rip)
    .globl cuGraphNodeGetDependentNodes_v2
    .type cuGraphNodeGetDependentNodes_v2,@function
cuGraphNodeGetDependentNodes_v2:
    jmp *kt_real_cuGraphNodeGetDependentNodes_v2(%rip)
    .globl cuGraphNodeGetEnabled
    .type cuGraphNodeGetEnabled,@function
cuGraphNodeGetEnabled:
    jmp *kt_real_cuGraphNodeGetEnabled(%rip)
    .globl cuGraphNodeGetLocalId
    .type cuGraphNodeGetLocalId,@function
cuGraphNodeGetLocalId:
    jmp *kt_real_cuGraphNodeGetLocalId(%rip)
    .globl cuGraphNodeGetParams
    .type cuGraphNodeGetParams,@function
cuGraphNodeGetParams:
    jmp *kt_real_cuGraphNodeGetParams(%rip)
    .globl cuGraphNodeGetToolsId
    .type cuGraphNodeGetToolsId,@function
cuGraphNodeGetToolsId:
    jmp *kt_real_cuGraphNodeGetToolsId(%rip)
    .globl cuGraphNodeGetType
    .type cuGraphNodeGetType,@function
cuGraphNodeGetType:
    jmp *kt_real_cuGraphNodeGetType(%rip)
    .globl cuGraphNodeSetEnabled
    .type cuGraphNodeSetEnabled,@function
cuGraphNodeSetEnabled:
    jmp *kt_real_cuGraphNodeSetEnabled(%rip)
    .globl cuGraphNodeSetParams
    .type cuGraphNodeSetParams,@function
cuGraphNodeSetParams:
    jmp *kt_real_cuGraphNodeSetParams(%rip)
    .globl cuGraphNodeSetParams_v2
    .type cuGraphNodeSetParams_v2,@function
cuGraphNodeSetParams_v2:
    jmp *kt_real_cuGraphNodeSetParams_v2(%rip)
    .globl cuGraphReleaseUserObject
    .type cuGraphReleaseUserObject,@function
cuGraphReleaseUserObject:
    jmp *kt_real_cuGraphReleaseUserObject(%rip)
    .globl cuGraphRemoveDependencies
    .type cuGraphRemoveDependencies,@function
cuGraphRemoveDependencies:
    jmp *kt_real_cuGraphRemoveDependencies(%rip)
    .globl cuGraphRemoveDependencies_v2
    .type cuGraphRemoveDependencies_v2,@function
cuGraphRemoveDependencies_v2:
    jmp *kt_real_cuGraphRemoveDependencies_v2(%rip)
    .globl cuGraphRetainUserObject
    .type cuGraphRetainUserObject,@function
cuGraphRetainUserObject:
    jmp *kt_real_cuGraphRetainUserObject(%rip)
    .globl cuGraphUpload
    .type cuGraphUpload,@function
cuGraphUpload:
    jmp *kt_real_cuGraphUpload(%rip)
    .globl cuGraphUpload_ptsz
    .type cuGraphUpload_ptsz,@function
cuGraphUpload_ptsz:
    jmp *kt_real_cuGraphUpload_ptsz(%rip)
    .globl cuGreenCtxCreate
    .type cuGreenCtxCreate,@function
cuGreenCtxCreate:
    jmp *kt_real_cuGreenCtxCreate(%rip)
    .globl cuGreenCtxDestroy
    .type cuGreenCtxDestroy,@function
cuGreenCtxDestroy:
    jmp *kt_real_cuGreenCtxDestroy(%rip)
    .globl cuGreenCtxGetDevResource
    .type cuGreenCtxGetDevResource,@function
cuGreenCtxGetDevResource:
    jmp *kt_real_cuGreenCtxGetDevResource(%rip)
    .globl cuGreenCtxGetId
    .type cuGreenCtxGetId,@function
cuGreenCtxGetId:
    jmp *kt_real_cuGreenCtxGetId(%rip)
    .globl cuGreenCtxRecordEvent
    .type cuGreenCtxRecordEvent,@function
cuGreenCtxRecordEvent:
    jmp *kt_real_cuGreenCtxRecordEvent(%rip)
    .globl cuGreenCtxStreamCreate
    .type cuGreenCtxStreamCreate,@function
cuGreenCtxStreamCreate:
    jmp *kt_real_cuGreenCtxStreamCreate(%rip)
    .globl cuGreenCtxWaitEvent
    .type cuGreenCtxWaitEvent,@function
cuGreenCtxWaitEvent:
    jmp *kt_real_cuGreenCtxWaitEvent(%rip)
    .globl cuImportExternalMemory
    .type cuImportExternalMemory,@function
cuImportExternalMemory:
    jmp *kt_real_cuImportExternalMemory(%rip)
    .globl cuImportExternalSemaphore
    .type cuImportExternalSemaphore,@function
cuImportExternalSemaphore:
    jmp *kt_real_cuImportExternalSemaphore(%rip)
    .globl cuInit
    .type cuInit,@function
cuInit:
    jmp *kt_real_cuInit(%rip)
    .globl cuIpcCloseMemHandle
    .type cuIpcCloseMemHandle,@function
cuIpcCloseMemHandle:
    jmp *kt_real_cuIpcCloseMemHandle(%rip)
    .globl cuIpcGetEventHandle
    .type cuIpcGetEventHandle,@function
cuIpcGetEventHandle:
    jmp *kt_real_cuIpcGetEventHandle(%rip)
    .globl cuIpcGetMemHandle
    .type cuIpcGetMemHandle,@function
cuIpcGetMemHandle:
    jmp *kt_real_cuIpcGetMemHandle(%rip)
    .globl cuIpcOpenEventHandle
    .type cuIpcOpenEventHandle,@function
cuIpcOpenEventHandle:
    jmp *kt_real_cuIpcOpenEventHandle(%rip)
    .globl cuIpcOpenMemHandle
    .type cuIpcOpenMemHandle,@function
cuIpcOpenMemHandle:
    jmp *kt_real_cuIpcOpenMemHandle(%rip)
    .globl cuIpcOpenMemHandle_v2
    .type cuIpcOpenMemHandle_v2,@function
cuIpcOpenMemHandle_v2:
    jmp *kt_real_cuIpcOpenMemHandle_v2(%rip)
    .globl cuKernelGetAttribute
    .type cuKernelGetAttribute,@function
cuKernelGetAttribute:
    jmp *kt_real_cuKernelGetAttribute(%rip)
    .globl cuKernelGetFunction
    .type cuKernelGetFunction,@function
cuKernelGetFunction:
    jmp *kt_real_cuKernelGetFunction(%rip)
    .globl cuKernelGetLibrary
    .type cuKernelGetLibrary,@function
cuKernelGetLibrary:
    jmp *kt_real_cuKernelGetLibrary(%rip)
    .globl cuKernelGetName
    .type cuKernelGetName,@function
cuKernelGetName:
    jmp *kt_real_cuKernelGetName(%rip)
    .globl cuKernelGetParamCount
    .type cuKernelGetParamCount,@function
cuKernelGetParamCount:
    jmp *kt_real_cuKernelGetParamCount(%rip)
    .globl cuKernelGetParamInfo
    .type cuKernelGetParamInfo,@function
cuKernelGetParamInfo:
    jmp *kt_real_cuKernelGetParamInfo(%rip)
    .globl cuKernelSetAttribute
    .type cuKernelSetAttribute,@function
cuKernelSetAttribute:
    jmp *kt_real_cuKernelSetAttribute(%rip)
    .globl cuKernelSetCacheConfig
    .type cuKernelSetCacheConfig,@function
cuKernelSetCacheConfig:
    jmp *kt_real_cuKernelSetCacheConfig(%rip)
    .globl cuLaunch
    .type cuLaunch,@function
cuLaunch:
    jmp *kt_real_cuLaunch(%rip)
    .globl cuLaunchCooperativeKernel
    .type cuLaunchCooperativeKernel,@function
cuLaunchCooperativeKernel:
    jmp *kt_real_cuLaunchCooperativeKernel(%rip)
    .globl cuLaunchCooperativeKernelMultiDevice
    .type cuLaunchCooperativeKernelMultiDevice,@function
cuLaunchCooperativeKernelMultiDevice:
    jmp *kt_real_cuLaunchCooperativeKernelMultiDevice(%rip)
    .globl cuLaunchCooperativeKernel_ptsz
    .type cuLaunchCooperativeKernel_ptsz,@function
cuLaunchCooperativeKernel_ptsz:
    jmp *kt_real_cuLaunchCooperativeKernel_ptsz(%rip)
    .globl cuLaunchGrid
    .type cuLaunchGrid,@function
cuLaunchGrid:
    jmp *kt_real_cuLaunchGrid(%rip)
    .globl cuLaunchGridAsync
    .type cuLaunchGridAsync,@function
cuLaunchGridAsync:
    jmp *kt_real_cuLaunchGridAsync(%rip)
    .globl cuLaunchHostFunc
    .type cuLaunchHostFunc,@function
cuLaunchHostFunc:
    jmp *kt_real_cuLaunchHostFunc(%rip)
    .globl cuLaunchHostFunc_ptsz
    .type cuLaunchHostFunc_ptsz,@function
cuLaunchHostFunc_ptsz:
    jmp *kt_real_cuLaunchHostFunc_ptsz(%rip)
    .globl cuLaunchHostFunc_v2
    .type cuLaunchHostFunc_v2,@function
cuLaunchHostFunc_v2:
    jmp *kt_real_cuLaunchHostFunc_v2(%rip)
    .globl cuLaunchHostFunc_v2_ptsz
    .type cuLaunchHostFunc_v2_ptsz,@function
cuLaunchHostFunc_v2_ptsz:
    jmp *kt_real_cuLaunchHostFunc_v2_ptsz(%rip)
    .globl cuLaunchKernelEx
    .type cuLaunchKernelEx,@function
cuLaunchKernelEx:
    jmp *kt_real_cuLaunchKernelEx(%rip)
    .globl cuLaunchKernelEx_ptsz
    .type cuLaunchKernelEx_ptsz,@function
cuLaunchKernelEx_ptsz:
    jmp *kt_real_cuLaunchKernelEx_ptsz(%rip)
    .globl cuLaunchKernel_ptsz
    .type cuLaunchKernel_ptsz,@function
cuLaunchKernel_ptsz:
    jmp *kt_real_cuLaunchKernel_ptsz(%rip)
    .globl cuLibraryEnumerateKernels
    .type cuLibraryEnumerateKernels,@function
cuLibraryEnumerateKernels:
    jmp *kt_real_cuLibraryEnumerateKernels(%rip)
    .globl cuLibraryGetGlobal
    .type cuLibraryGetGlobal,@function
cuLibraryGetGlobal:
    jmp *kt_real_cuLibraryGetGlobal(%rip)
    .globl cuLibraryGetKernel
    .type cuLibraryGetKernel,@function
cuLibraryGetKernel:
    jmp *kt_real_cuLibraryGetKernel(%rip)
    .globl cuLibraryGetKernelCount
    .type cuLibraryGetKernelCount,@function
cuLibraryGetKernelCount:
    jmp *kt_real_cuLibraryGetKernelCount(%rip)
    .globl cuLibraryGetManaged
    .type cuLibraryGetManaged,@function
cuLibraryGetManaged:
    jmp *kt_real_cuLibraryGetManaged(%rip)
    .globl cuLibraryGetModule
    .type cuLibraryGetModule,@function
cuLibraryGetModule:
    jmp *kt_real_cuLibraryGetModule(%rip)
    .globl cuLibraryGetUnifiedFunction
    .type cuLibraryGetUnifiedFunction,@function
cuLibraryGetUnifiedFunction:
    jmp *kt_real_cuLibraryGetUnifiedFunction(%rip)
    .globl cuLibraryLoadData
    .type cuLibraryLoadData,@function
cuLibraryLoadData:
    jmp *kt_real_cuLibraryLoadData(%rip)
    .globl cuLibraryLoadFromFile
    .type cuLibraryLoadFromFile,@function
cuLibraryLoadFromFile:
    jmp *kt_real_cuLibraryLoadFromFile(%rip)
    .globl cuLibraryUnload
    .type cuLibraryUnload,@function
cuLibraryUnload:
    jmp *kt_real_cuLibraryUnload(%rip)
    .globl cuLinkAddData
    .type cuLinkAddData,@function
cuLinkAddData:
    jmp *kt_real_cuLinkAddData(%rip)
    .globl cuLinkAddData_v2
    .type cuLinkAddData_v2,@function
cuLinkAddData_v2:
    jmp *kt_real_cuLinkAddData_v2(%rip)
    .globl cuLinkAddFile
    .type cuLinkAddFile,@function
cuLinkAddFile:
    jmp *kt_real_cuLinkAddFile(%rip)
    .globl cuLinkAddFile_v2
    .type cuLinkAddFile_v2,@function
cuLinkAddFile_v2:
    jmp *kt_real_cuLinkAddFile_v2(%rip)
    .globl cuLinkComplete
    .type cuLinkComplete,@function
cuLinkComplete:
    jmp *kt_real_cuLinkComplete(%rip)
    .globl cuLinkCreate
    .type cuLinkCreate,@function
cuLinkCreate:
    jmp *kt_real_cuLinkCreate(%rip)
    .globl cuLinkCreate_v2
    .type cuLinkCreate_v2,@function
cuLinkCreate_v2:
    jmp *kt_real_cuLinkCreate_v2(%rip)
    .globl cuLinkDestroy
    .type cuLinkDestroy,@function
cuLinkDestroy:
    jmp *kt_real_cuLinkDestroy(%rip)
    .globl cuLogicalEndpointAddDevice
    .type cuLogicalEndpointAddDevice,@function
cuLogicalEndpointAddDevice:
    jmp *kt_real_cuLogicalEndpointAddDevice(%rip)
    .globl cuLogicalEndpointBindAddr
    .type cuLogicalEndpointBindAddr,@function
cuLogicalEndpointBindAddr:
    jmp *kt_real_cuLogicalEndpointBindAddr(%rip)
    .globl cuLogicalEndpointBindMem
    .type cuLogicalEndpointBindMem,@function
cuLogicalEndpointBindMem:
    jmp *kt_real_cuLogicalEndpointBindMem(%rip)
    .globl cuLogicalEndpointCreate
    .type cuLogicalEndpointCreate,@function
cuLogicalEndpointCreate:
    jmp *kt_real_cuLogicalEndpointCreate(%rip)
    .globl cuLogicalEndpointDestroy
    .type cuLogicalEndpointDestroy,@function
cuLogicalEndpointDestroy:
    jmp *kt_real_cuLogicalEndpointDestroy(%rip)
    .globl cuLogicalEndpointExport
    .type cuLogicalEndpointExport,@function
cuLogicalEndpointExport:
    jmp *kt_real_cuLogicalEndpointExport(%rip)
    .globl cuLogicalEndpointGetLimits
    .type cuLogicalEndpointGetLimits,@function
cuLogicalEndpointGetLimits:
    jmp *kt_real_cuLogicalEndpointGetLimits(%rip)
    .globl cuLogicalEndpointIdRelease
    .type cuLogicalEndpointIdRelease,@function
cuLogicalEndpointIdRelease:
    jmp *kt_real_cuLogicalEndpointIdRelease(%rip)
    .globl cuLogicalEndpointIdReserve
    .type cuLogicalEndpointIdReserve,@function
cuLogicalEndpointIdReserve:
    jmp *kt_real_cuLogicalEndpointIdReserve(%rip)
    .globl cuLogicalEndpointImport
    .type cuLogicalEndpointImport,@function
cuLogicalEndpointImport:
    jmp *kt_real_cuLogicalEndpointImport(%rip)
    .globl cuLogicalEndpointQuery
    .type cuLogicalEndpointQuery,@function
cuLogicalEndpointQuery:
    jmp *kt_real_cuLogicalEndpointQuery(%rip)
    .globl cuLogicalEndpointUnbind
    .type cuLogicalEndpointUnbind,@function
cuLogicalEndpointUnbind:
    jmp *kt_real_cuLogicalEndpointUnbind(%rip)
    .globl cuLogsCurrent
    .type cuLogsCurrent,@function
cuLogsCurrent:
    jmp *kt_real_cuLogsCurrent(%rip)
    .globl cuLogsDumpToFile
    .type cuLogsDumpToFile,@function
cuLogsDumpToFile:
    jmp *kt_real_cuLogsDumpToFile(%rip)
    .globl cuLogsDumpToMemory
    .type cuLogsDumpToMemory,@function
cuLogsDumpToMemory:
    jmp *kt_real_cuLogsDumpToMemory(%rip)
    .globl cuLogsRegisterCallback
    .type cuLogsRegisterCallback,@function
cuLogsRegisterCallback:
    jmp *kt_real_cuLogsRegisterCallback(%rip)
    .globl cuLogsUnregisterCallback
    .type cuLogsUnregisterCallback,@function
cuLogsUnregisterCallback:
    jmp *kt_real_cuLogsUnregisterCallback(%rip)
    .globl cuMemAddressFree
    .type cuMemAddressFree,@function
cuMemAddressFree:
    jmp *kt_real_cuMemAddressFree(%rip)
    .globl cuMemAddressReserve
    .type cuMemAddressReserve,@function
cuMemAddressReserve:
    jmp *kt_real_cuMemAddressReserve(%rip)
    .globl cuMemAdvise
    .type cuMemAdvise,@function
cuMemAdvise:
    jmp *kt_real_cuMemAdvise(%rip)
    .globl cuMemAdvise_v2
    .type cuMemAdvise_v2,@function
cuMemAdvise_v2:
    jmp *kt_real_cuMemAdvise_v2(%rip)
    .globl cuMemAlloc
    .type cuMemAlloc,@function
cuMemAlloc:
    jmp *kt_real_cuMemAlloc(%rip)
    .globl cuMemAllocAsync
    .type cuMemAllocAsync,@function
cuMemAllocAsync:
    jmp *kt_real_cuMemAllocAsync(%rip)
    .globl cuMemAllocAsync_ptsz
    .type cuMemAllocAsync_ptsz,@function
cuMemAllocAsync_ptsz:
    jmp *kt_real_cuMemAllocAsync_ptsz(%rip)
    .globl cuMemAllocFromPoolAsync
    .type cuMemAllocFromPoolAsync,@function
cuMemAllocFromPoolAsync:
    jmp *kt_real_cuMemAllocFromPoolAsync(%rip)
    .globl cuMemAllocFromPoolAsync_ptsz
    .type cuMemAllocFromPoolAsync_ptsz,@function
cuMemAllocFromPoolAsync_ptsz:
    jmp *kt_real_cuMemAllocFromPoolAsync_ptsz(%rip)
    .globl cuMemAllocHost
    .type cuMemAllocHost,@function
cuMemAllocHost:
    jmp *kt_real_cuMemAllocHost(%rip)
    .globl cuMemAllocHost_v2
    .type cuMemAllocHost_v2,@function
cuMemAllocHost_v2:
    jmp *kt_real_cuMemAllocHost_v2(%rip)
    .globl cuMemAllocManaged
    .type cuMemAllocManaged,@function
cuMemAllocManaged:
    jmp *kt_real_cuMemAllocManaged(%rip)
    .globl cuMemAllocPitch
    .type cuMemAllocPitch,@function
cuMemAllocPitch:
    jmp *kt_real_cuMemAllocPitch(%rip)
    .globl cuMemAllocPitch_v2
    .type cuMemAllocPitch_v2,@function
cuMemAllocPitch_v2:
    jmp *kt_real_cuMemAllocPitch_v2(%rip)
    .globl cuMemAlloc_v2
    .type cuMemAlloc_v2,@function
cuMemAlloc_v2:
    jmp *kt_real_cuMemAlloc_v2(%rip)
    .globl cuMemBatchDecompressAsync
    .type cuMemBatchDecompressAsync,@function
cuMemBatchDecompressAsync:
    jmp *kt_real_cuMemBatchDecompressAsync(%rip)
    .globl cuMemBatchDecompressAsync_ptsz
    .type cuMemBatchDecompressAsync_ptsz,@function
cuMemBatchDecompressAsync_ptsz:
    jmp *kt_real_cuMemBatchDecompressAsync_ptsz(%rip)
    .globl cuMemcpy
    .type cuMemcpy,@function
cuMemcpy:
    jmp *kt_real_cuMemcpy(%rip)
    .globl cuMemcpy2D
    .type cuMemcpy2D,@function
cuMemcpy2D:
    jmp *kt_real_cuMemcpy2D(%rip)
    .globl cuMemcpy2DAsync
    .type cuMemcpy2DAsync,@function
cuMemcpy2DAsync:
    jmp *kt_real_cuMemcpy2DAsync(%rip)
    .globl cuMemcpy2DAsync_v2
    .type cuMemcpy2DAsync_v2,@function
cuMemcpy2DAsync_v2:
    jmp *kt_real_cuMemcpy2DAsync_v2(%rip)
    .globl cuMemcpy2DAsync_v2_ptsz
    .type cuMemcpy2DAsync_v2_ptsz,@function
cuMemcpy2DAsync_v2_ptsz:
    jmp *kt_real_cuMemcpy2DAsync_v2_ptsz(%rip)
    .globl cuMemcpy2DUnaligned
    .type cuMemcpy2DUnaligned,@function
cuMemcpy2DUnaligned:
    jmp *kt_real_cuMemcpy2DUnaligned(%rip)
    .globl cuMemcpy2DUnaligned_v2
    .type cuMemcpy2DUnaligned_v2,@function
cuMemcpy2DUnaligned_v2:
    jmp *kt_real_cuMemcpy2DUnaligned_v2(%rip)
    .globl cuMemcpy2DUnaligned_v2_ptds
    .type cuMemcpy2DUnaligned_v2_ptds,@function
cuMemcpy2DUnaligned_v2_ptds:
    jmp *kt_real_cuMemcpy2DUnaligned_v2_ptds(%rip)
    .globl cuMemcpy2D_v2
    .type cuMemcpy2D_v2,@function
cuMemcpy2D_v2:
    jmp *kt_real_cuMemcpy2D_v2(%rip)
    .globl cuMemcpy2D_v2_ptds
    .type cuMemcpy2D_v2_ptds,@function
cuMemcpy2D_v2_ptds:
    jmp *kt_real_cuMemcpy2D_v2_ptds(%rip)
    .globl cuMemcpy3D
    .type cuMemcpy3D,@function
cuMemcpy3D:
    jmp *kt_real_cuMemcpy3D(%rip)
    .globl cuMemcpy3DAsync
    .type cuMemcpy3DAsync,@function
cuMemcpy3DAsync:
    jmp *kt_real_cuMemcpy3DAsync(%rip)
    .globl cuMemcpy3DAsync_v2
    .type cuMemcpy3DAsync_v2,@function
cuMemcpy3DAsync_v2:
    jmp *kt_real_cuMemcpy3DAsync_v2(%rip)
    .globl cuMemcpy3DAsync_v2_ptsz
    .type cuMemcpy3DAsync_v2_ptsz,@function
cuMemcpy3DAsync_v2_ptsz:
    jmp *kt_real_cuMemcpy3DAsync_v2_ptsz(%rip)
    .globl cuMemcpy3DBatchAsync
    .type cuMemcpy3DBatchAsync,@function
cuMemcpy3DBatchAsync:
    jmp *kt_real_cuMemcpy3DBatchAsync(%rip)
    .globl cuMemcpy3DBatchAsync_ptsz
    .type cuMemcpy3DBatchAsync_ptsz,@function
cuMemcpy3DBatchAsync_ptsz:
    jmp *kt_real_cuMemcpy3DBatchAsync_ptsz(%rip)
    .globl cuMemcpy3DBatchAsync_v2
    .type cuMemcpy3DBatchAsync_v2,@function
cuMemcpy3DBatchAsync_v2:
    jmp *kt_real_cuMemcpy3DBatchAsync_v2(%rip)
    .globl cuMemcpy3DBatchAsync_v2_ptsz
    .type cuMemcpy3DBatchAsync_v2_ptsz,@function
cuMemcpy3DBatchAsync_v2_ptsz:
    jmp *kt_real_cuMemcpy3DBatchAsync_v2_ptsz(%rip)
    .globl cuMemcpy3DPeer
    .type cuMemcpy3DPeer,@function
cuMemcpy3DPeer:
    jmp *kt_real_cuMemcpy3DPeer(%rip)
    .globl cuMemcpy3DPeerAsync
    .type cuMemcpy3DPeerAsync,@function
cuMemcpy3DPeerAsync:
    jmp *kt_real_cuMemcpy3DPeerAsync(%rip)
    .globl cuMemcpy3DPeerAsync_ptsz
    .type cuMemcpy3DPeerAsync_ptsz,@function
cuMemcpy3DPeerAsync_ptsz:
    jmp *kt_real_cuMemcpy3DPeerAsync_ptsz(%rip)
    .globl cuMemcpy3DPeer_ptds
    .type cuMemcpy3DPeer_ptds,@function
cuMemcpy3DPeer_ptds:
    jmp *kt_real_cuMemcpy3DPeer_ptds(%rip)
    .globl cuMemcpy3D_v2
    .type cuMemcpy3D_v2,@function
cuMemcpy3D_v2:
    jmp *kt_real_cuMemcpy3D_v2(%rip)
    .globl cuMemcpy3D_v2_ptds
    .type cuMemcpy3D_v2_ptds,@function
cuMemcpy3D_v2_ptds:
    jmp *kt_real_cuMemcpy3D_v2_ptds(%rip)
    .globl cuMemcpy3DWithAttributesAsync
    .type cuMemcpy3DWithAttributesAsync,@function
cuMemcpy3DWithAttributesAsync:
    jmp *kt_real_cuMemcpy3DWithAttributesAsync(%rip)
    .globl cuMemcpy3DWithAttributesAsync_ptsz
    .type cuMemcpy3DWithAttributesAsync_ptsz,@function
cuMemcpy3DWithAttributesAsync_ptsz:
    jmp *kt_real_cuMemcpy3DWithAttributesAsync_ptsz(%rip)
    .globl cuMemcpyAsync
    .type cuMemcpyAsync,@function
cuMemcpyAsync:
    jmp *kt_real_cuMemcpyAsync(%rip)
    .globl cuMemcpyAsync_ptsz
    .type cuMemcpyAsync_ptsz,@function
cuMemcpyAsync_ptsz:
    jmp *kt_real_cuMemcpyAsync_ptsz(%rip)
    .globl cuMemcpyAtoA
    .type cuMemcpyAtoA,@function
cuMemcpyAtoA:
    jmp *kt_real_cuMemcpyAtoA(%rip)
    .globl cuMemcpyAtoA_v2
    .type cuMemcpyAtoA_v2,@function
cuMemcpyAtoA_v2:
    jmp *kt_real_cuMemcpyAtoA_v2(%rip)
    .globl cuMemcpyAtoA_v2_ptds
    .type cuMemcpyAtoA_v2_ptds,@function
cuMemcpyAtoA_v2_ptds:
    jmp *kt_real_cuMemcpyAtoA_v2_ptds(%rip)
    .globl cuMemcpyAtoD
    .type cuMemcpyAtoD,@function
cuMemcpyAtoD:
    jmp *kt_real_cuMemcpyAtoD(%rip)
    .globl cuMemcpyAtoD_v2
    .type cuMemcpyAtoD_v2,@function
cuMemcpyAtoD_v2:
    jmp *kt_real_cuMemcpyAtoD_v2(%rip)
    .globl cuMemcpyAtoD_v2_ptds
    .type cuMemcpyAtoD_v2_ptds,@function
cuMemcpyAtoD_v2_ptds:
    jmp *kt_real_cuMemcpyAtoD_v2_ptds(%rip)
    .globl cuMemcpyAtoH
    .type cuMemcpyAtoH,@function
cuMemcpyAtoH:
    jmp *kt_real_cuMemcpyAtoH(%rip)
    .globl cuMemcpyAtoHAsync
    .type cuMemcpyAtoHAsync,@function
cuMemcpyAtoHAsync:
    jmp *kt_real_cuMemcpyAtoHAsync(%rip)
    .globl cuMemcpyAtoHAsync_v2
    .type cuMemcpyAtoHAsync_v2,@function
cuMemcpyAtoHAsync_v2:
    jmp *kt_real_cuMemcpyAtoHAsync_v2(%rip)
    .globl cuMemcpyAtoHAsync_v2_ptsz
    .type cuMemcpyAtoHAsync_v2_ptsz,@function
cuMemcpyAtoHAsync_v2_ptsz:
    jmp *kt_real_cuMemcpyAtoHAsync_v2_ptsz(%rip)
    .globl cuMemcpyAtoH_v2
    .type cuMemcpyAtoH_v2,@function
cuMemcpyAtoH_v2:
    jmp *kt_real_cuMemcpyAtoH_v2(%rip)
    .globl cuMemcpyAtoH_v2_ptds
    .type cuMemcpyAtoH_v2_ptds,@function
cuMemcpyAtoH_v2_ptds:
    jmp *kt_real_cuMemcpyAtoH_v2_ptds(%rip)
    .globl cuMemcpyBatchAsync
    .type cuMemcpyBatchAsync,@function
cuMemcpyBatchAsync:
    jmp *kt_real_cuMemcpyBatchAsync(%rip)
    .globl cuMemcpyBatchAsync_ptsz
    .type cuMemcpyBatchAsync_ptsz,@function
cuMemcpyBatchAsync_ptsz:
    jmp *kt_real_cuMemcpyBatchAsync_ptsz(%rip)
    .globl cuMemcpyBatchAsync_v2
    .type cuMemcpyBatchAsync_v2,@function
cuMemcpyBatchAsync_v2:
    jmp *kt_real_cuMemcpyBatchAsync_v2(%rip)
    .globl cuMemcpyBatchAsync_v2_ptsz
    .type cuMemcpyBatchAsync_v2_ptsz,@function
cuMemcpyBatchAsync_v2_ptsz:
    jmp *kt_real_cuMemcpyBatchAsync_v2_ptsz(%rip)
    .globl cuMemcpyDtoA
    .type cuMemcpyDtoA,@function
cuMemcpyDtoA:
    jmp *kt_real_cuMemcpyDtoA(%rip)
    .globl cuMemcpyDtoA_v2
    .type cuMemcpyDtoA_v2,@function
cuMemcpyDtoA_v2:
    jmp *kt_real_cuMemcpyDtoA_v2(%rip)
    .globl cuMemcpyDtoA_v2_ptds
    .type cuMemcpyDtoA_v2_ptds,@function
cuMemcpyDtoA_v2_ptds:
    jmp *kt_real_cuMemcpyDtoA_v2_ptds(%rip)
    .globl cuMemcpyDtoD
    .type cuMemcpyDtoD,@function
cuMemcpyDtoD:
    jmp *kt_real_cuMemcpyDtoD(%rip)
    .globl cuMemcpyDtoDAsync
    .type cuMemcpyDtoDAsync,@function
cuMemcpyDtoDAsync:
    jmp *kt_real_cuMemcpyDtoDAsync(%rip)
    .globl cuMemcpyDtoDAsync_v2
    .type cuMemcpyDtoDAsync_v2,@function
cuMemcpyDtoDAsync_v2:
    jmp *kt_real_cuMemcpyDtoDAsync_v2(%rip)
    .globl cuMemcpyDtoDAsync_v2_ptsz
    .type cuMemcpyDtoDAsync_v2_ptsz,@function
cuMemcpyDtoDAsync_v2_ptsz:
    jmp *kt_real_cuMemcpyDtoDAsync_v2_ptsz(%rip)
    .globl cuMemcpyDtoD_v2
    .type cuMemcpyDtoD_v2,@function
cuMemcpyDtoD_v2:
    jmp *kt_real_cuMemcpyDtoD_v2(%rip)
    .globl cuMemcpyDtoD_v2_ptds
    .type cuMemcpyDtoD_v2_ptds,@function
cuMemcpyDtoD_v2_ptds:
    jmp *kt_real_cuMemcpyDtoD_v2_ptds(%rip)
    .globl cuMemcpyDtoH
    .type cuMemcpyDtoH,@function
cuMemcpyDtoH:
    jmp *kt_real_cuMemcpyDtoH(%rip)
    .globl cuMemcpyDtoHAsync
    .type cuMemcpyDtoHAsync,@function
cuMemcpyDtoHAsync:
    jmp *kt_real_cuMemcpyDtoHAsync(%rip)
    .globl cuMemcpyDtoHAsync_v2
    .type cuMemcpyDtoHAsync_v2,@function
cuMemcpyDtoHAsync_v2:
    jmp *kt_real_cuMemcpyDtoHAsync_v2(%rip)
    .globl cuMemcpyDtoHAsync_v2_ptsz
    .type cuMemcpyDtoHAsync_v2_ptsz,@function
cuMemcpyDtoHAsync_v2_ptsz:
    jmp *kt_real_cuMemcpyDtoHAsync_v2_ptsz(%rip)
    .globl cuMemcpyDtoH_v2
    .type cuMemcpyDtoH_v2,@function
cuMemcpyDtoH_v2:
    jmp *kt_real_cuMemcpyDtoH_v2(%rip)
    .globl cuMemcpyDtoH_v2_ptds
    .type cuMemcpyDtoH_v2_ptds,@function
cuMemcpyDtoH_v2_ptds:
    jmp *kt_real_cuMemcpyDtoH_v2_ptds(%rip)
    .globl cuMemcpyHtoA
    .type cuMemcpyHtoA,@function
cuMemcpyHtoA:
    jmp *kt_real_cuMemcpyHtoA(%rip)
    .globl cuMemcpyHtoAAsync
    .type cuMemcpyHtoAAsync,@function
cuMemcpyHtoAAsync:
    jmp *kt_real_cuMemcpyHtoAAsync(%rip)
    .globl cuMemcpyHtoAAsync_v2
    .type cuMemcpyHtoAAsync_v2,@function
cuMemcpyHtoAAsync_v2:
    jmp *kt_real_cuMemcpyHtoAAsync_v2(%rip)
    .globl cuMemcpyHtoAAsync_v2_ptsz
    .type cuMemcpyHtoAAsync_v2_ptsz,@function
cuMemcpyHtoAAsync_v2_ptsz:
    jmp *kt_real_cuMemcpyHtoAAsync_v2_ptsz(%rip)
    .globl cuMemcpyHtoA_v2
    .type cuMemcpyHtoA_v2,@function
cuMemcpyHtoA_v2:
    jmp *kt_real_cuMemcpyHtoA_v2(%rip)
    .globl cuMemcpyHtoA_v2_ptds
    .type cuMemcpyHtoA_v2_ptds,@function
cuMemcpyHtoA_v2_ptds:
    jmp *kt_real_cuMemcpyHtoA_v2_ptds(%rip)
    .globl cuMemcpyHtoD
    .type cuMemcpyHtoD,@function
cuMemcpyHtoD:
    jmp *kt_real_cuMemcpyHtoD(%rip)
    .globl cuMemcpyHtoDAsync
    .type cuMemcpyHtoDAsync,@function
cuMemcpyHtoDAsync:
    jmp *kt_real_cuMemcpyHtoDAsync(%rip)
    .globl cuMemcpyHtoDAsync_v2
    .type cuMemcpyHtoDAsync_v2,@function
cuMemcpyHtoDAsync_v2:
    jmp *kt_real_cuMemcpyHtoDAsync_v2(%rip)
    .globl cuMemcpyHtoDAsync_v2_ptsz
    .type cuMemcpyHtoDAsync_v2_ptsz,@function
cuMemcpyHtoDAsync_v2_ptsz:
    jmp *kt_real_cuMemcpyHtoDAsync_v2_ptsz(%rip)
    .globl cuMemcpyHtoD_v2
    .type cuMemcpyHtoD_v2,@function
cuMemcpyHtoD_v2:
    jmp *kt_real_cuMemcpyHtoD_v2(%rip)
    .globl cuMemcpyHtoD_v2_ptds
    .type cuMemcpyHtoD_v2_ptds,@function
cuMemcpyHtoD_v2_ptds:
    jmp *kt_real_cuMemcpyHtoD_v2_ptds(%rip)
    .globl cuMemcpyPeer
    .type cuMemcpyPeer,@function
cuMemcpyPeer:
    jmp *kt_real_cuMemcpyPeer(%rip)
    .globl cuMemcpyPeerAsync
    .type cuMemcpyPeerAsync,@function
cuMemcpyPeerAsync:
    jmp *kt_real_cuMemcpyPeerAsync(%rip)
    .globl cuMemcpyPeerAsync_ptsz
    .type cuMemcpyPeerAsync_ptsz,@function
cuMemcpyPeerAsync_ptsz:
    jmp *kt_real_cuMemcpyPeerAsync_ptsz(%rip)
    .globl cuMemcpyPeer_ptds
    .type cuMemcpyPeer_ptds,@function
cuMemcpyPeer_ptds:
    jmp *kt_real_cuMemcpyPeer_ptds(%rip)
    .globl cuMemcpy_ptds
    .type cuMemcpy_ptds,@function
cuMemcpy_ptds:
    jmp *kt_real_cuMemcpy_ptds(%rip)
    .globl cuMemcpyWithAttributesAsync
    .type cuMemcpyWithAttributesAsync,@function
cuMemcpyWithAttributesAsync:
    jmp *kt_real_cuMemcpyWithAttributesAsync(%rip)
    .globl cuMemcpyWithAttributesAsync_ptsz
    .type cuMemcpyWithAttributesAsync_ptsz,@function
cuMemcpyWithAttributesAsync_ptsz:
    jmp *kt_real_cuMemcpyWithAttributesAsync_ptsz(%rip)
    .globl cuMemCreate
    .type cuMemCreate,@function
cuMemCreate:
    jmp *kt_real_cuMemCreate(%rip)
    .globl cuMemDiscardAndPrefetchBatchAsync
    .type cuMemDiscardAndPrefetchBatchAsync,@function
cuMemDiscardAndPrefetchBatchAsync:
    jmp *kt_real_cuMemDiscardAndPrefetchBatchAsync(%rip)
    .globl cuMemDiscardAndPrefetchBatchAsync_ptsz
    .type cuMemDiscardAndPrefetchBatchAsync_ptsz,@function
cuMemDiscardAndPrefetchBatchAsync_ptsz:
    jmp *kt_real_cuMemDiscardAndPrefetchBatchAsync_ptsz(%rip)
    .globl cuMemDiscardBatchAsync
    .type cuMemDiscardBatchAsync,@function
cuMemDiscardBatchAsync:
    jmp *kt_real_cuMemDiscardBatchAsync(%rip)
    .globl cuMemDiscardBatchAsync_ptsz
    .type cuMemDiscardBatchAsync_ptsz,@function
cuMemDiscardBatchAsync_ptsz:
    jmp *kt_real_cuMemDiscardBatchAsync_ptsz(%rip)
    .globl cuMemExportToShareableHandle
    .type cuMemExportToShareableHandle,@function
cuMemExportToShareableHandle:
    jmp *kt_real_cuMemExportToShareableHandle(%rip)
    .globl cuMemFree
    .type cuMemFree,@function
cuMemFree:
    jmp *kt_real_cuMemFree(%rip)
    .globl cuMemFreeAsync
    .type cuMemFreeAsync,@function
cuMemFreeAsync:
    jmp *kt_real_cuMemFreeAsync(%rip)
    .globl cuMemFreeAsync_ptsz
    .type cuMemFreeAsync_ptsz,@function
cuMemFreeAsync_ptsz:
    jmp *kt_real_cuMemFreeAsync_ptsz(%rip)
    .globl cuMemFreeHost
    .type cuMemFreeHost,@function
cuMemFreeHost:
    jmp *kt_real_cuMemFreeHost(%rip)
    .globl cuMemFree_v2
    .type cuMemFree_v2,@function
cuMemFree_v2:
    jmp *kt_real_cuMemFree_v2(%rip)
    .globl cuMemGetAccess
    .type cuMemGetAccess,@function
cuMemGetAccess:
    jmp *kt_real_cuMemGetAccess(%rip)
    .globl cuMemGetAddressRange
    .type cuMemGetAddressRange,@function
cuMemGetAddressRange:
    jmp *kt_real_cuMemGetAddressRange(%rip)
    .globl cuMemGetAddressRange_v2
    .type cuMemGetAddressRange_v2,@function
cuMemGetAddressRange_v2:
    jmp *kt_real_cuMemGetAddressRange_v2(%rip)
    .globl cuMemGetAllocationGranularity
    .type cuMemGetAllocationGranularity,@function
cuMemGetAllocationGranularity:
    jmp *kt_real_cuMemGetAllocationGranularity(%rip)
    .globl cuMemGetAllocationPropertiesFromHandle
    .type cuMemGetAllocationPropertiesFromHandle,@function
cuMemGetAllocationPropertiesFromHandle:
    jmp *kt_real_cuMemGetAllocationPropertiesFromHandle(%rip)
    .globl cuMemGetAttribute
    .type cuMemGetAttribute,@function
cuMemGetAttribute:
    jmp *kt_real_cuMemGetAttribute(%rip)
    .globl cuMemGetAttribute_v2
    .type cuMemGetAttribute_v2,@function
cuMemGetAttribute_v2:
    jmp *kt_real_cuMemGetAttribute_v2(%rip)
    .globl cuMemGetDefaultMemPool
    .type cuMemGetDefaultMemPool,@function
cuMemGetDefaultMemPool:
    jmp *kt_real_cuMemGetDefaultMemPool(%rip)
    .globl cuMemGetHandleForAddressRange
    .type cuMemGetHandleForAddressRange,@function
cuMemGetHandleForAddressRange:
    jmp *kt_real_cuMemGetHandleForAddressRange(%rip)
    .globl cuMemGetInfo
    .type cuMemGetInfo,@function
cuMemGetInfo:
    jmp *kt_real_cuMemGetInfo(%rip)
    .globl cuMemGetInfo_v2
    .type cuMemGetInfo_v2,@function
cuMemGetInfo_v2:
    jmp *kt_real_cuMemGetInfo_v2(%rip)
    .globl cuMemGetLocationInfo
    .type cuMemGetLocationInfo,@function
cuMemGetLocationInfo:
    jmp *kt_real_cuMemGetLocationInfo(%rip)
    .globl cuMemGetMemPool
    .type cuMemGetMemPool,@function
cuMemGetMemPool:
    jmp *kt_real_cuMemGetMemPool(%rip)
    .globl cuMemHostAlloc
    .type cuMemHostAlloc,@function
cuMemHostAlloc:
    jmp *kt_real_cuMemHostAlloc(%rip)
    .globl cuMemHostGetDevicePointer
    .type cuMemHostGetDevicePointer,@function
cuMemHostGetDevicePointer:
    jmp *kt_real_cuMemHostGetDevicePointer(%rip)
    .globl cuMemHostGetDevicePointer_v2
    .type cuMemHostGetDevicePointer_v2,@function
cuMemHostGetDevicePointer_v2:
    jmp *kt_real_cuMemHostGetDevicePointer_v2(%rip)
    .globl cuMemHostGetFlags
    .type cuMemHostGetFlags,@function
cuMemHostGetFlags:
    jmp *kt_real_cuMemHostGetFlags(%rip)
    .globl cuMemHostRegister
    .type cuMemHostRegister,@function
cuMemHostRegister:
    jmp *kt_real_cuMemHostRegister(%rip)
    .globl cuMemHostRegister_v2
    .type cuMemHostRegister_v2,@function
cuMemHostRegister_v2:
    jmp *kt_real_cuMemHostRegister_v2(%rip)
    .globl cuMemHostUnregister
    .type cuMemHostUnregister,@function
cuMemHostUnregister:
    jmp *kt_real_cuMemHostUnregister(%rip)
    .globl cuMemImportFromShareableHandle
    .type cuMemImportFromShareableHandle,@function
cuMemImportFromShareableHandle:
    jmp *kt_real_cuMemImportFromShareableHandle(%rip)
    .globl cuMemMap
    .type cuMemMap,@function
cuMemMap:
    jmp *kt_real_cuMemMap(%rip)
    .globl cuMemMapArrayAsync
    .type cuMemMapArrayAsync,@function
cuMemMapArrayAsync:
    jmp *kt_real_cuMemMapArrayAsync(%rip)
    .globl cuMemMapArrayAsync_ptsz
    .type cuMemMapArrayAsync_ptsz,@function
cuMemMapArrayAsync_ptsz:
    jmp *kt_real_cuMemMapArrayAsync_ptsz(%rip)
    .globl cuMemPoolCreate
    .type cuMemPoolCreate,@function
cuMemPoolCreate:
    jmp *kt_real_cuMemPoolCreate(%rip)
    .globl cuMemPoolDestroy
    .type cuMemPoolDestroy,@function
cuMemPoolDestroy:
    jmp *kt_real_cuMemPoolDestroy(%rip)
    .globl cuMemPoolExportPointer
    .type cuMemPoolExportPointer,@function
cuMemPoolExportPointer:
    jmp *kt_real_cuMemPoolExportPointer(%rip)
    .globl cuMemPoolExportToShareableHandle
    .type cuMemPoolExportToShareableHandle,@function
cuMemPoolExportToShareableHandle:
    jmp *kt_real_cuMemPoolExportToShareableHandle(%rip)
    .globl cuMemPoolGetAccess
    .type cuMemPoolGetAccess,@function
cuMemPoolGetAccess:
    jmp *kt_real_cuMemPoolGetAccess(%rip)
    .globl cuMemPoolGetAttribute
    .type cuMemPoolGetAttribute,@function
cuMemPoolGetAttribute:
    jmp *kt_real_cuMemPoolGetAttribute(%rip)
    .globl cuMemPoolImportFromShareableHandle
    .type cuMemPoolImportFromShareableHandle,@function
cuMemPoolImportFromShareableHandle:
    jmp *kt_real_cuMemPoolImportFromShareableHandle(%rip)
    .globl cuMemPoolImportPointer
    .type cuMemPoolImportPointer,@function
cuMemPoolImportPointer:
    jmp *kt_real_cuMemPoolImportPointer(%rip)
    .globl cuMemPoolSetAccess
    .type cuMemPoolSetAccess,@function
cuMemPoolSetAccess:
    jmp *kt_real_cuMemPoolSetAccess(%rip)
    .globl cuMemPoolSetAttribute
    .type cuMemPoolSetAttribute,@function
cuMemPoolSetAttribute:
    jmp *kt_real_cuMemPoolSetAttribute(%rip)
    .globl cuMemPoolTrimTo
    .type cuMemPoolTrimTo,@function
cuMemPoolTrimTo:
    jmp *kt_real_cuMemPoolTrimTo(%rip)
    .globl cuMemPrefetchAsync
    .type cuMemPrefetchAsync,@function
cuMemPrefetchAsync:
    jmp *kt_real_cuMemPrefetchAsync(%rip)
    .globl cuMemPrefetchAsync_ptsz
    .type cuMemPrefetchAsync_ptsz,@function
cuMemPrefetchAsync_ptsz:
    jmp *kt_real_cuMemPrefetchAsync_ptsz(%rip)
    .globl cuMemPrefetchAsync_v2
    .type cuMemPrefetchAsync_v2,@function
cuMemPrefetchAsync_v2:
    jmp *kt_real_cuMemPrefetchAsync_v2(%rip)
    .globl cuMemPrefetchAsync_v2_ptsz
    .type cuMemPrefetchAsync_v2_ptsz,@function
cuMemPrefetchAsync_v2_ptsz:
    jmp *kt_real_cuMemPrefetchAsync_v2_ptsz(%rip)
    .globl cuMemPrefetchBatchAsync
    .type cuMemPrefetchBatchAsync,@function
cuMemPrefetchBatchAsync:
    jmp *kt_real_cuMemPrefetchBatchAsync(%rip)
    .globl cuMemPrefetchBatchAsync_ptsz
    .type cuMemPrefetchBatchAsync_ptsz,@function
cuMemPrefetchBatchAsync_ptsz:
    jmp *kt_real_cuMemPrefetchBatchAsync_ptsz(%rip)
    .globl cuMemRangeGetAttribute
    .type cuMemRangeGetAttribute,@function
cuMemRangeGetAttribute:
    jmp *kt_real_cuMemRangeGetAttribute(%rip)
    .globl cuMemRangeGetAttributes
    .type cuMemRangeGetAttributes,@function
cuMemRangeGetAttributes:
    jmp *kt_real_cuMemRangeGetAttributes(%rip)
    .globl cuMemRelease
    .type cuMemRelease,@function
cuMemRelease:
    jmp *kt_real_cuMemRelease(%rip)
    .globl cuMemRetainAllocationHandle
    .type cuMemRetainAllocationHandle,@function
cuMemRetainAllocationHandle:
    jmp *kt_real_cuMemRetainAllocationHandle(%rip)
    .globl cuMemSetAccess
    .type cuMemSetAccess,@function
cuMemSetAccess:
    jmp *kt_real_cuMemSetAccess(%rip)
    .globl cuMemsetD16
    .type cuMemsetD16,@function
cuMemsetD16:
    jmp *kt_real_cuMemsetD16(%rip)
    .globl cuMemsetD16Async
    .type cuMemsetD16Async,@function
cuMemsetD16Async:
    jmp *kt_real_cuMemsetD16Async(%rip)
    .globl cuMemsetD16Async_ptsz
    .type cuMemsetD16Async_ptsz,@function
cuMemsetD16Async_ptsz:
    jmp *kt_real_cuMemsetD16Async_ptsz(%rip)
    .globl cuMemsetD16_v2
    .type cuMemsetD16_v2,@function
cuMemsetD16_v2:
    jmp *kt_real_cuMemsetD16_v2(%rip)
    .globl cuMemsetD16_v2_ptds
    .type cuMemsetD16_v2_ptds,@function
cuMemsetD16_v2_ptds:
    jmp *kt_real_cuMemsetD16_v2_ptds(%rip)
    .globl cuMemsetD2D16
    .type cuMemsetD2D16,@function
cuMemsetD2D16:
    jmp *kt_real_cuMemsetD2D16(%rip)
    .globl cuMemsetD2D16Async
    .type cuMemsetD2D16Async,@function
cuMemsetD2D16Async:
    jmp *kt_real_cuMemsetD2D16Async(%rip)
    .globl cuMemsetD2D16Async_ptsz
    .type cuMemsetD2D16Async_ptsz,@function
cuMemsetD2D16Async_ptsz:
    jmp *kt_real_cuMemsetD2D16Async_ptsz(%rip)
    .globl cuMemsetD2D16_v2
    .type cuMemsetD2D16_v2,@function
cuMemsetD2D16_v2:
    jmp *kt_real_cuMemsetD2D16_v2(%rip)
    .globl cuMemsetD2D16_v2_ptds
    .type cuMemsetD2D16_v2_ptds,@function
cuMemsetD2D16_v2_ptds:
    jmp *kt_real_cuMemsetD2D16_v2_ptds(%rip)
    .globl cuMemsetD2D32
    .type cuMemsetD2D32,@function
cuMemsetD2D32:
    jmp *kt_real_cuMemsetD2D32(%rip)
    .globl cuMemsetD2D32Async
    .type cuMemsetD2D32Async,@function
cuMemsetD2D32Async:
    jmp *kt_real_cuMemsetD2D32Async(%rip)
    .globl cuMemsetD2D32Async_ptsz
    .type cuMemsetD2D32Async_ptsz,@function
cuMemsetD2D32Async_ptsz:
    jmp *kt_real_cuMemsetD2D32Async_ptsz(%rip)
    .globl cuMemsetD2D32_v2
    .type cuMemsetD2D32_v2,@function
cuMemsetD2D32_v2:
    jmp *kt_real_cuMemsetD2D32_v2(%rip)
    .globl cuMemsetD2D32_v2_ptds
    .type cuMemsetD2D32_v2_ptds,@function
cuMemsetD2D32_v2_ptds:
    jmp *kt_real_cuMemsetD2D32_v2_ptds(%rip)
    .globl cuMemsetD2D8
    .type cuMemsetD2D8,@function
cuMemsetD2D8:
    jmp *kt_real_cuMemsetD2D8(%rip)
    .globl cuMemsetD2D8Async
    .type cuMemsetD2D8Async,@function
cuMemsetD2D8Async:
    jmp *kt_real_cuMemsetD2D8Async(%rip)
    .globl cuMemsetD2D8Async_ptsz
    .type cuMemsetD2D8Async_ptsz,@function
cuMemsetD2D8Async_ptsz:
    jmp *kt_real_cuMemsetD2D8Async_ptsz(%rip)
    .globl cuMemsetD2D8_v2
    .type cuMemsetD2D8_v2,@function
cuMemsetD2D8_v2:
    jmp *kt_real_cuMemsetD2D8_v2(%rip)
    .globl cuMemsetD2D8_v2_ptds
    .type cuMemsetD2D8_v2_ptds,@function
cuMemsetD2D8_v2_ptds:
    jmp *kt_real_cuMemsetD2D8_v2_ptds(%rip)
    .globl cuMemsetD32
    .type cuMemsetD32,@function
cuMemsetD32:
    jmp *kt_real_cuMemsetD32(%rip)
    .globl cuMemsetD32Async
    .type cuMemsetD32Async,@function
cuMemsetD32Async:
    jmp *kt_real_cuMemsetD32Async(%rip)
    .globl cuMemsetD32Async_ptsz
    .type cuMemsetD32Async_ptsz,@function
cuMemsetD32Async_ptsz:
    jmp *kt_real_cuMemsetD32Async_ptsz(%rip)
    .globl cuMemsetD32_v2
    .type cuMemsetD32_v2,@function
cuMemsetD32_v2:
    jmp *kt_real_cuMemsetD32_v2(%rip)
    .globl cuMemsetD32_v2_ptds
    .type cuMemsetD32_v2_ptds,@function
cuMemsetD32_v2_ptds:
    jmp *kt_real_cuMemsetD32_v2_ptds(%rip)
    .globl cuMemsetD8
    .type cuMemsetD8,@function
cuMemsetD8:
    jmp *kt_real_cuMemsetD8(%rip)
    .globl cuMemsetD8Async
    .type cuMemsetD8Async,@function
cuMemsetD8Async:
    jmp *kt_real_cuMemsetD8Async(%rip)
    .globl cuMemsetD8Async_ptsz
    .type cuMemsetD8Async_ptsz,@function
cuMemsetD8Async_ptsz:
    jmp *kt_real_cuMemsetD8Async_ptsz(%rip)
    .globl cuMemsetD8_v2
    .type cuMemsetD8_v2,@function
cuMemsetD8_v2:
    jmp *kt_real_cuMemsetD8_v2(%rip)
    .globl cuMemsetD8_v2_ptds
    .type cuMemsetD8_v2_ptds,@function
cuMemsetD8_v2_ptds:
    jmp *kt_real_cuMemsetD8_v2_ptds(%rip)
    .globl cuMemSetMemPool
    .type cuMemSetMemPool,@function
cuMemSetMemPool:
    jmp *kt_real_cuMemSetMemPool(%rip)
    .globl cuMemUnmap
    .type cuMemUnmap,@function
cuMemUnmap:
    jmp *kt_real_cuMemUnmap(%rip)
    .globl cuMipmappedArrayCreate
    .type cuMipmappedArrayCreate,@function
cuMipmappedArrayCreate:
    jmp *kt_real_cuMipmappedArrayCreate(%rip)
    .globl cuMipmappedArrayDestroy
    .type cuMipmappedArrayDestroy,@function
cuMipmappedArrayDestroy:
    jmp *kt_real_cuMipmappedArrayDestroy(%rip)
    .globl cuMipmappedArrayGetLevel
    .type cuMipmappedArrayGetLevel,@function
cuMipmappedArrayGetLevel:
    jmp *kt_real_cuMipmappedArrayGetLevel(%rip)
    .globl cuMipmappedArrayGetMemoryRequirements
    .type cuMipmappedArrayGetMemoryRequirements,@function
cuMipmappedArrayGetMemoryRequirements:
    jmp *kt_real_cuMipmappedArrayGetMemoryRequirements(%rip)
    .globl cuMipmappedArrayGetSparseProperties
    .type cuMipmappedArrayGetSparseProperties,@function
cuMipmappedArrayGetSparseProperties:
    jmp *kt_real_cuMipmappedArrayGetSparseProperties(%rip)
    .globl cuModuleEnumerateFunctions
    .type cuModuleEnumerateFunctions,@function
cuModuleEnumerateFunctions:
    jmp *kt_real_cuModuleEnumerateFunctions(%rip)
    .globl cuModuleGetFunction
    .type cuModuleGetFunction,@function
cuModuleGetFunction:
    jmp *kt_real_cuModuleGetFunction(%rip)
    .globl cuModuleGetFunctionCount
    .type cuModuleGetFunctionCount,@function
cuModuleGetFunctionCount:
    jmp *kt_real_cuModuleGetFunctionCount(%rip)
    .globl cuModuleGetGlobal
    .type cuModuleGetGlobal,@function
cuModuleGetGlobal:
    jmp *kt_real_cuModuleGetGlobal(%rip)
    .globl cuModuleGetGlobal_v2
    .type cuModuleGetGlobal_v2,@function
cuModuleGetGlobal_v2:
    jmp *kt_real_cuModuleGetGlobal_v2(%rip)
    .globl cuModuleGetLoadingMode
    .type cuModuleGetLoadingMode,@function
cuModuleGetLoadingMode:
    jmp *kt_real_cuModuleGetLoadingMode(%rip)
    .globl cuModuleGetSurfRef
    .type cuModuleGetSurfRef,@function
cuModuleGetSurfRef:
    jmp *kt_real_cuModuleGetSurfRef(%rip)
    .globl cuModuleGetTexRef
    .type cuModuleGetTexRef,@function
cuModuleGetTexRef:
    jmp *kt_real_cuModuleGetTexRef(%rip)
    .globl cuModuleLoad
    .type cuModuleLoad,@function
cuModuleLoad:
    jmp *kt_real_cuModuleLoad(%rip)
    .globl cuModuleLoadData
    .type cuModuleLoadData,@function
cuModuleLoadData:
    jmp *kt_real_cuModuleLoadData(%rip)
    .globl cuModuleLoadDataEx
    .type cuModuleLoadDataEx,@function
cuModuleLoadDataEx:
    jmp *kt_real_cuModuleLoadDataEx(%rip)
    .globl cuModuleLoadFatBinary
    .type cuModuleLoadFatBinary,@function
cuModuleLoadFatBinary:
    jmp *kt_real_cuModuleLoadFatBinary(%rip)
    .globl cuModuleUnload
    .type cuModuleUnload,@function
cuModuleUnload:
    jmp *kt_real_cuModuleUnload(%rip)
    .globl cuMulticastAddDevice
    .type cuMulticastAddDevice,@function
cuMulticastAddDevice:
    jmp *kt_real_cuMulticastAddDevice(%rip)
    .globl cuMulticastBindAddr
    .type cuMulticastBindAddr,@function
cuMulticastBindAddr:
    jmp *kt_real_cuMulticastBindAddr(%rip)
    .globl cuMulticastBindAddr_v2
    .type cuMulticastBindAddr_v2,@function
cuMulticastBindAddr_v2:
    jmp *kt_real_cuMulticastBindAddr_v2(%rip)
    .globl cuMulticastBindMem
    .type cuMulticastBindMem,@function
cuMulticastBindMem:
    jmp *kt_real_cuMulticastBindMem(%rip)
    .globl cuMulticastBindMem_v2
    .type cuMulticastBindMem_v2,@function
cuMulticastBindMem_v2:
    jmp *kt_real_cuMulticastBindMem_v2(%rip)
    .globl cuMulticastCreate
    .type cuMulticastCreate,@function
cuMulticastCreate:
    jmp *kt_real_cuMulticastCreate(%rip)
    .globl cuMulticastGetGranularity
    .type cuMulticastGetGranularity,@function
cuMulticastGetGranularity:
    jmp *kt_real_cuMulticastGetGranularity(%rip)
    .globl cuMulticastUnbind
    .type cuMulticastUnbind,@function
cuMulticastUnbind:
    jmp *kt_real_cuMulticastUnbind(%rip)
    .globl cuOccupancyAvailableDynamicSMemPerBlock
    .type cuOccupancyAvailableDynamicSMemPerBlock,@function
cuOccupancyAvailableDynamicSMemPerBlock:
    jmp *kt_real_cuOccupancyAvailableDynamicSMemPerBlock(%rip)
    .globl cuOccupancyMaxActiveBlocksPerMultiprocessor
    .type cuOccupancyMaxActiveBlocksPerMultiprocessor,@function
cuOccupancyMaxActiveBlocksPerMultiprocessor:
    jmp *kt_real_cuOccupancyMaxActiveBlocksPerMultiprocessor(%rip)
    .globl cuOccupancyMaxActiveBlocksPerMultiprocessorWithFlags
    .type cuOccupancyMaxActiveBlocksPerMultiprocessorWithFlags,@function
cuOccupancyMaxActiveBlocksPerMultiprocessorWithFlags:
    jmp *kt_real_cuOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(%rip)
    .globl cuOccupancyMaxActiveClusters
    .type cuOccupancyMaxActiveClusters,@function
cuOccupancyMaxActiveClusters:
    jmp *kt_real_cuOccupancyMaxActiveClusters(%rip)
    .globl cuOccupancyMaxPotentialBlockSize
    .type cuOccupancyMaxPotentialBlockSize,@function
cuOccupancyMaxPotentialBlockSize:
    jmp *kt_real_cuOccupancyMaxPotentialBlockSize(%rip)
    .globl cuOccupancyMaxPotentialBlockSizeWithFlags
    .type cuOccupancyMaxPotentialBlockSizeWithFlags,@function
cuOccupancyMaxPotentialBlockSizeWithFlags:
    jmp *kt_real_cuOccupancyMaxPotentialBlockSizeWithFlags(%rip)
    .globl cuOccupancyMaxPotentialClusterSize
    .type cuOccupancyMaxPotentialClusterSize,@function
cuOccupancyMaxPotentialClusterSize:
    jmp *kt_real_cuOccupancyMaxPotentialClusterSize(%rip)
    .globl cuParamSetf
    .type cuParamSetf,@function
cuParamSetf:
    jmp *kt_real_cuParamSetf(%rip)
    .globl cuParamSeti
    .type cuParamSeti,@function
cuParamSeti:
    jmp *kt_real_cuParamSeti(%rip)
    .globl cuParamSetSize
    .type cuParamSetSize,@function
cuParamSetSize:
    jmp *kt_real_cuParamSetSize(%rip)
    .globl cuParamSetTexRef
    .type cuParamSetTexRef,@function
cuParamSetTexRef:
    jmp *kt_real_cuParamSetTexRef(%rip)
    .globl cuParamSetv
    .type cuParamSetv,@function
cuParamSetv:
    jmp *kt_real_cuParamSetv(%rip)
    .globl cuPointerGetAttribute
    .type cuPointerGetAttribute,@function
cuPointerGetAttribute:
    jmp *kt_real_cuPointerGetAttribute(%rip)
    .globl cuPointerGetAttributes
    .type cuPointerGetAttributes,@function
cuPointerGetAttributes:
    jmp *kt_real_cuPointerGetAttributes(%rip)
    .globl cuPointerSetAttribute
    .type cuPointerSetAttribute,@function
cuPointerSetAttribute:
    jmp *kt_real_cuPointerSetAttribute(%rip)
    .globl cuProfilerInitialize
    .type cuProfilerInitialize,@function
cuProfilerInitialize:
    jmp *kt_real_cuProfilerInitialize(%rip)
    .globl cuProfilerStart
    .type cuProfilerStart,@function
cuProfilerStart:
    jmp *kt_real_cuProfilerStart(%rip)
    .globl cuProfilerStop
    .type cuProfilerStop,@function
cuProfilerStop:
    jmp *kt_real_cuProfilerStop(%rip)
    .globl cuSignalExternalSemaphoresAsync
    .type cuSignalExternalSemaphoresAsync,@function
cuSignalExternalSemaphoresAsync:
    jmp *kt_real_cuSignalExternalSemaphoresAsync(%rip)
    .globl cuSignalExternalSemaphoresAsync_ptsz
    .type cuSignalExternalSemaphoresAsync_ptsz,@function
cuSignalExternalSemaphoresAsync_ptsz:
    jmp *kt_real_cuSignalExternalSemaphoresAsync_ptsz(%rip)
    .globl cuStreamAddCallback
    .type cuStreamAddCallback,@function
cuStreamAddCallback:
    jmp *kt_real_cuStreamAddCallback(%rip)
    .globl cuStreamAddCallback_ptsz
    .type cuStreamAddCallback_ptsz,@function
cuStreamAddCallback_ptsz:
    jmp *kt_real_cuStreamAddCallback_ptsz(%rip)
    .globl cuStreamAttachMemAsync
    .type cuStreamAttachMemAsync,@function
cuStreamAttachMemAsync:
    jmp *kt_real_cuStreamAttachMemAsync(%rip)
    .globl cuStreamAttachMemAsync_ptsz
    .type cuStreamAttachMemAsync_ptsz,@function
cuStreamAttachMemAsync_ptsz:
    jmp *kt_real_cuStreamAttachMemAsync_ptsz(%rip)
    .globl cuStreamBatchMemOp
    .type cuStreamBatchMemOp,@function
cuStreamBatchMemOp:
    jmp *kt_real_cuStreamBatchMemOp(%rip)
    .globl cuStreamBatchMemOp_ptsz
    .type cuStreamBatchMemOp_ptsz,@function
cuStreamBatchMemOp_ptsz:
    jmp *kt_real_cuStreamBatchMemOp_ptsz(%rip)
    .globl cuStreamBatchMemOp_v2
    .type cuStreamBatchMemOp_v2,@function
cuStreamBatchMemOp_v2:
    jmp *kt_real_cuStreamBatchMemOp_v2(%rip)
    .globl cuStreamBatchMemOp_v2_ptsz
    .type cuStreamBatchMemOp_v2_ptsz,@function
cuStreamBatchMemOp_v2_ptsz:
    jmp *kt_real_cuStreamBatchMemOp_v2_ptsz(%rip)
    .globl cuStreamBeginCapture
    .type cuStreamBeginCapture,@function
cuStreamBeginCapture:
    jmp *kt_real_cuStreamBeginCapture(%rip)
    .globl cuStreamBeginCapture_ptsz
    .type cuStreamBeginCapture_ptsz,@function
cuStreamBeginCapture_ptsz:
    jmp *kt_real_cuStreamBeginCapture_ptsz(%rip)
    .globl cuStreamBeginCaptureToCig
    .type cuStreamBeginCaptureToCig,@function
cuStreamBeginCaptureToCig:
    jmp *kt_real_cuStreamBeginCaptureToCig(%rip)
    .globl cuStreamBeginCaptureToCig_ptsz
    .type cuStreamBeginCaptureToCig_ptsz,@function
cuStreamBeginCaptureToCig_ptsz:
    jmp *kt_real_cuStreamBeginCaptureToCig_ptsz(%rip)
    .globl cuStreamBeginCaptureToGraph
    .type cuStreamBeginCaptureToGraph,@function
cuStreamBeginCaptureToGraph:
    jmp *kt_real_cuStreamBeginCaptureToGraph(%rip)
    .globl cuStreamBeginCaptureToGraph_ptsz
    .type cuStreamBeginCaptureToGraph_ptsz,@function
cuStreamBeginCaptureToGraph_ptsz:
    jmp *kt_real_cuStreamBeginCaptureToGraph_ptsz(%rip)
    .globl cuStreamBeginCapture_v2
    .type cuStreamBeginCapture_v2,@function
cuStreamBeginCapture_v2:
    jmp *kt_real_cuStreamBeginCapture_v2(%rip)
    .globl cuStreamBeginCapture_v2_ptsz
    .type cuStreamBeginCapture_v2_ptsz,@function
cuStreamBeginCapture_v2_ptsz:
    jmp *kt_real_cuStreamBeginCapture_v2_ptsz(%rip)
    .globl cuStreamBeginRecaptureToGraph
    .type cuStreamBeginRecaptureToGraph,@function
cuStreamBeginRecaptureToGraph:
    jmp *kt_real_cuStreamBeginRecaptureToGraph(%rip)
    .globl cuStreamBeginRecaptureToGraph_ptsz
    .type cuStreamBeginRecaptureToGraph_ptsz,@function
cuStreamBeginRecaptureToGraph_ptsz:
    jmp *kt_real_cuStreamBeginRecaptureToGraph_ptsz(%rip)
    .globl cuStreamCopyAttributes
    .type cuStreamCopyAttributes,@function
cuStreamCopyAttributes:
    jmp *kt_real_cuStreamCopyAttributes(%rip)
    .globl cuStreamCopyAttributes_ptsz
    .type cuStreamCopyAttributes_ptsz,@function
cuStreamCopyAttributes_ptsz:
    jmp *kt_real_cuStreamCopyAttributes_ptsz(%rip)
    .globl cuStreamCreate
    .type cuStreamCreate,@function
cuStreamCreate:
    jmp *kt_real_cuStreamCreate(%rip)
    .globl cuStreamCreateWithPriority
    .type cuStreamCreateWithPriority,@function
cuStreamCreateWithPriority:
    jmp *kt_real_cuStreamCreateWithPriority(%rip)
    .globl cuStreamDestroy
    .type cuStreamDestroy,@function
cuStreamDestroy:
    jmp *kt_real_cuStreamDestroy(%rip)
    .globl cuStreamDestroy_v2
    .type cuStreamDestroy_v2,@function
cuStreamDestroy_v2:
    jmp *kt_real_cuStreamDestroy_v2(%rip)
    .globl cuStreamEndCapture
    .type cuStreamEndCapture,@function
cuStreamEndCapture:
    jmp *kt_real_cuStreamEndCapture(%rip)
    .globl cuStreamEndCapture_ptsz
    .type cuStreamEndCapture_ptsz,@function
cuStreamEndCapture_ptsz:
    jmp *kt_real_cuStreamEndCapture_ptsz(%rip)
    .globl cuStreamEndCaptureToCig
    .type cuStreamEndCaptureToCig,@function
cuStreamEndCaptureToCig:
    jmp *kt_real_cuStreamEndCaptureToCig(%rip)
    .globl cuStreamEndCaptureToCig_ptsz
    .type cuStreamEndCaptureToCig_ptsz,@function
cuStreamEndCaptureToCig_ptsz:
    jmp *kt_real_cuStreamEndCaptureToCig_ptsz(%rip)
    .globl cuStreamGetAttribute
    .type cuStreamGetAttribute,@function
cuStreamGetAttribute:
    jmp *kt_real_cuStreamGetAttribute(%rip)
    .globl cuStreamGetAttribute_ptsz
    .type cuStreamGetAttribute_ptsz,@function
cuStreamGetAttribute_ptsz:
    jmp *kt_real_cuStreamGetAttribute_ptsz(%rip)
    .globl cuStreamGetCaptureInfo
    .type cuStreamGetCaptureInfo,@function
cuStreamGetCaptureInfo:
    jmp *kt_real_cuStreamGetCaptureInfo(%rip)
    .globl cuStreamGetCaptureInfo_ptsz
    .type cuStreamGetCaptureInfo_ptsz,@function
cuStreamGetCaptureInfo_ptsz:
    jmp *kt_real_cuStreamGetCaptureInfo_ptsz(%rip)
    .globl cuStreamGetCaptureInfo_v2
    .type cuStreamGetCaptureInfo_v2,@function
cuStreamGetCaptureInfo_v2:
    jmp *kt_real_cuStreamGetCaptureInfo_v2(%rip)
    .globl cuStreamGetCaptureInfo_v2_ptsz
    .type cuStreamGetCaptureInfo_v2_ptsz,@function
cuStreamGetCaptureInfo_v2_ptsz:
    jmp *kt_real_cuStreamGetCaptureInfo_v2_ptsz(%rip)
    .globl cuStreamGetCaptureInfo_v3
    .type cuStreamGetCaptureInfo_v3,@function
cuStreamGetCaptureInfo_v3:
    jmp *kt_real_cuStreamGetCaptureInfo_v3(%rip)
    .globl cuStreamGetCaptureInfo_v3_ptsz
    .type cuStreamGetCaptureInfo_v3_ptsz,@function
cuStreamGetCaptureInfo_v3_ptsz:
    jmp *kt_real_cuStreamGetCaptureInfo_v3_ptsz(%rip)
    .globl cuStreamGetCtx
    .type cuStreamGetCtx,@function
cuStreamGetCtx:
    jmp *kt_real_cuStreamGetCtx(%rip)
    .globl cuStreamGetCtx_ptsz
    .type cuStreamGetCtx_ptsz,@function
cuStreamGetCtx_ptsz:
    jmp *kt_real_cuStreamGetCtx_ptsz(%rip)
    .globl cuStreamGetCtx_v2
    .type cuStreamGetCtx_v2,@function
cuStreamGetCtx_v2:
    jmp *kt_real_cuStreamGetCtx_v2(%rip)
    .globl cuStreamGetCtx_v2_ptsz
    .type cuStreamGetCtx_v2_ptsz,@function
cuStreamGetCtx_v2_ptsz:
    jmp *kt_real_cuStreamGetCtx_v2_ptsz(%rip)
    .globl cuStreamGetDevice
    .type cuStreamGetDevice,@function
cuStreamGetDevice:
    jmp *kt_real_cuStreamGetDevice(%rip)
    .globl cuStreamGetDevice_ptsz
    .type cuStreamGetDevice_ptsz,@function
cuStreamGetDevice_ptsz:
    jmp *kt_real_cuStreamGetDevice_ptsz(%rip)
    .globl cuStreamGetDevResource
    .type cuStreamGetDevResource,@function
cuStreamGetDevResource:
    jmp *kt_real_cuStreamGetDevResource(%rip)
    .globl cuStreamGetDevResource_ptsz
    .type cuStreamGetDevResource_ptsz,@function
cuStreamGetDevResource_ptsz:
    jmp *kt_real_cuStreamGetDevResource_ptsz(%rip)
    .globl cuStreamGetFlags
    .type cuStreamGetFlags,@function
cuStreamGetFlags:
    jmp *kt_real_cuStreamGetFlags(%rip)
    .globl cuStreamGetFlags_ptsz
    .type cuStreamGetFlags_ptsz,@function
cuStreamGetFlags_ptsz:
    jmp *kt_real_cuStreamGetFlags_ptsz(%rip)
    .globl cuStreamGetGreenCtx
    .type cuStreamGetGreenCtx,@function
cuStreamGetGreenCtx:
    jmp *kt_real_cuStreamGetGreenCtx(%rip)
    .globl cuStreamGetId
    .type cuStreamGetId,@function
cuStreamGetId:
    jmp *kt_real_cuStreamGetId(%rip)
    .globl cuStreamGetId_ptsz
    .type cuStreamGetId_ptsz,@function
cuStreamGetId_ptsz:
    jmp *kt_real_cuStreamGetId_ptsz(%rip)
    .globl cuStreamGetPriority
    .type cuStreamGetPriority,@function
cuStreamGetPriority:
    jmp *kt_real_cuStreamGetPriority(%rip)
    .globl cuStreamGetPriority_ptsz
    .type cuStreamGetPriority_ptsz,@function
cuStreamGetPriority_ptsz:
    jmp *kt_real_cuStreamGetPriority_ptsz(%rip)
    .globl cuStreamIsCapturing
    .type cuStreamIsCapturing,@function
cuStreamIsCapturing:
    jmp *kt_real_cuStreamIsCapturing(%rip)
    .globl cuStreamIsCapturing_ptsz
    .type cuStreamIsCapturing_ptsz,@function
cuStreamIsCapturing_ptsz:
    jmp *kt_real_cuStreamIsCapturing_ptsz(%rip)
    .globl cuStreamQuery
    .type cuStreamQuery,@function
cuStreamQuery:
    jmp *kt_real_cuStreamQuery(%rip)
    .globl cuStreamQuery_ptsz
    .type cuStreamQuery_ptsz,@function
cuStreamQuery_ptsz:
    jmp *kt_real_cuStreamQuery_ptsz(%rip)
    .globl cuStreamSetAttribute
    .type cuStreamSetAttribute,@function
cuStreamSetAttribute:
    jmp *kt_real_cuStreamSetAttribute(%rip)
    .globl cuStreamSetAttribute_ptsz
    .type cuStreamSetAttribute_ptsz,@function
cuStreamSetAttribute_ptsz:
    jmp *kt_real_cuStreamSetAttribute_ptsz(%rip)
    .globl cuStreamSynchronize_ptsz
    .type cuStreamSynchronize_ptsz,@function
cuStreamSynchronize_ptsz:
    jmp *kt_real_cuStreamSynchronize_ptsz(%rip)
    .globl cuStreamUpdateCaptureDependencies
    .type cuStreamUpdateCaptureDependencies,@function
cuStreamUpdateCaptureDependencies:
    jmp *kt_real_cuStreamUpdateCaptureDependencies(%rip)
    .globl cuStreamUpdateCaptureDependencies_ptsz
    .type cuStreamUpdateCaptureDependencies_ptsz,@function
cuStreamUpdateCaptureDependencies_ptsz:
    jmp *kt_real_cuStreamUpdateCaptureDependencies_ptsz(%rip)
    .globl cuStreamUpdateCaptureDependencies_v2
    .type cuStreamUpdateCaptureDependencies_v2,@function
cuStreamUpdateCaptureDependencies_v2:
    jmp *kt_real_cuStreamUpdateCaptureDependencies_v2(%rip)
    .globl cuStreamUpdateCaptureDependencies_v2_ptsz
    .type cuStreamUpdateCaptureDependencies_v2_ptsz,@function
cuStreamUpdateCaptureDependencies_v2_ptsz:
    jmp *kt_real_cuStreamUpdateCaptureDependencies_v2_ptsz(%rip)
    .globl cuStreamWaitEvent
    .type cuStreamWaitEvent,@function
cuStreamWaitEvent:
    jmp *kt_real_cuStreamWaitEvent(%rip)
    .globl cuStreamWaitEvent_ptsz
    .type cuStreamWaitEvent_ptsz,@function
cuStreamWaitEvent_ptsz:
    jmp *kt_real_cuStreamWaitEvent_ptsz(%rip)
    .globl cuStreamWaitValue32
    .type cuStreamWaitValue32,@function
cuStreamWaitValue32:
    jmp *kt_real_cuStreamWaitValue32(%rip)
    .globl cuStreamWaitValue32_ptsz
    .type cuStreamWaitValue32_ptsz,@function
cuStreamWaitValue32_ptsz:
    jmp *kt_real_cuStreamWaitValue32_ptsz(%rip)
    .globl cuStreamWaitValue32_v2
    .type cuStreamWaitValue32_v2,@function
cuStreamWaitValue32_v2:
    jmp *kt_real_cuStreamWaitValue32_v2(%rip)
    .globl cuStreamWaitValue32_v2_ptsz
    .type cuStreamWaitValue32_v2_ptsz,@function
cuStreamWaitValue32_v2_ptsz:
    jmp *kt_real_cuStreamWaitValue32_v2_ptsz(%rip)
    .globl cuStreamWaitValue64
    .type cuStreamWaitValue64,@function
cuStreamWaitValue64:
    jmp *kt_real_cuStreamWaitValue64(%rip)
    .globl cuStreamWaitValue64_ptsz
    .type cuStreamWaitValue64_ptsz,@function
cuStreamWaitValue64_ptsz:
    jmp *kt_real_cuStreamWaitValue64_ptsz(%rip)
    .globl cuStreamWaitValue64_v2
    .type cuStreamWaitValue64_v2,@function
cuStreamWaitValue64_v2:
    jmp *kt_real_cuStreamWaitValue64_v2(%rip)
    .globl cuStreamWaitValue64_v2_ptsz
    .type cuStreamWaitValue64_v2_ptsz,@function
cuStreamWaitValue64_v2_ptsz:
    jmp *kt_real_cuStreamWaitValue64_v2_ptsz(%rip)
    .globl cuStreamWriteValue32
    .type cuStreamWriteValue32,@function
cuStreamWriteValue32:
    jmp *kt_real_cuStreamWriteValue32(%rip)
    .globl cuStreamWriteValue32_ptsz
    .type cuStreamWriteValue32_ptsz,@function
cuStreamWriteValue32_ptsz:
    jmp *kt_real_cuStreamWriteValue32_ptsz(%rip)
    .globl cuStreamWriteValue32_v2
    .type cuStreamWriteValue32_v2,@function
cuStreamWriteValue32_v2:
    jmp *kt_real_cuStreamWriteValue32_v2(%rip)
    .globl cuStreamWriteValue32_v2_ptsz
    .type cuStreamWriteValue32_v2_ptsz,@function
cuStreamWriteValue32_v2_ptsz:
    jmp *kt_real_cuStreamWriteValue32_v2_ptsz(%rip)
    .globl cuStreamWriteValue64
    .type cuStreamWriteValue64,@function
cuStreamWriteValue64:
    jmp *kt_real_cuStreamWriteValue64(%rip)
    .globl cuStreamWriteValue64_ptsz
    .type cuStreamWriteValue64_ptsz,@function
cuStreamWriteValue64_ptsz:
    jmp *kt_real_cuStreamWriteValue64_ptsz(%rip)
    .globl cuStreamWriteValue64_v2
    .type cuStreamWriteValue64_v2,@function
cuStreamWriteValue64_v2:
    jmp *kt_real_cuStreamWriteValue64_v2(%rip)
    .globl cuStreamWriteValue64_v2_ptsz
    .type cuStreamWriteValue64_v2_ptsz,@function
cuStreamWriteValue64_v2_ptsz:
    jmp *kt_real_cuStreamWriteValue64_v2_ptsz(%rip)
    .globl cuSubgridCreate
    .type cuSubgridCreate,@function
cuSubgridCreate:
    jmp *kt_real_cuSubgridCreate(%rip)
    .globl cuSubgridWorkerGridCreate
    .type cuSubgridWorkerGridCreate,@function
cuSubgridWorkerGridCreate:
    jmp *kt_real_cuSubgridWorkerGridCreate(%rip)
    .globl cuSubgridWorksetCreate
    .type cuSubgridWorksetCreate,@function
cuSubgridWorksetCreate:
    jmp *kt_real_cuSubgridWorksetCreate(%rip)
    .globl cuSurfObjectCreate
    .type cuSurfObjectCreate,@function
cuSurfObjectCreate:
    jmp *kt_real_cuSurfObjectCreate(%rip)
    .globl cuSurfObjectDestroy
    .type cuSurfObjectDestroy,@function
cuSurfObjectDestroy:
    jmp *kt_real_cuSurfObjectDestroy(%rip)
    .globl cuSurfObjectGetResourceDesc
    .type cuSurfObjectGetResourceDesc,@function
cuSurfObjectGetResourceDesc:
    jmp *kt_real_cuSurfObjectGetResourceDesc(%rip)
    .globl cuSurfRefGetArray
    .type cuSurfRefGetArray,@function
cuSurfRefGetArray:
    jmp *kt_real_cuSurfRefGetArray(%rip)
    .globl cuSurfRefSetArray
    .type cuSurfRefSetArray,@function
cuSurfRefSetArray:
    jmp *kt_real_cuSurfRefSetArray(%rip)
    .globl cuTensorMapEncodeIm2col
    .type cuTensorMapEncodeIm2col,@function
cuTensorMapEncodeIm2col:
    jmp *kt_real_cuTensorMapEncodeIm2col(%rip)
    .globl cuTensorMapEncodeIm2colWide
    .type cuTensorMapEncodeIm2colWide,@function
cuTensorMapEncodeIm2colWide:
    jmp *kt_real_cuTensorMapEncodeIm2colWide(%rip)
    .globl cuTensorMapEncodeTiled
    .type cuTensorMapEncodeTiled,@function
cuTensorMapEncodeTiled:
    jmp *kt_real_cuTensorMapEncodeTiled(%rip)
    .globl cuTensorMapReplaceAddress
    .type cuTensorMapReplaceAddress,@function
cuTensorMapReplaceAddress:
    jmp *kt_real_cuTensorMapReplaceAddress(%rip)
    .globl cuTexObjectCreate
    .type cuTexObjectCreate,@function
cuTexObjectCreate:
    jmp *kt_real_cuTexObjectCreate(%rip)
    .globl cuTexObjectDestroy
    .type cuTexObjectDestroy,@function
cuTexObjectDestroy:
    jmp *kt_real_cuTexObjectDestroy(%rip)
    .globl cuTexObjectGetResourceDesc
    .type cuTexObjectGetResourceDesc,@function
cuTexObjectGetResourceDesc:
    jmp *kt_real_cuTexObjectGetResourceDesc(%rip)
    .globl cuTexObjectGetResourceViewDesc
    .type cuTexObjectGetResourceViewDesc,@function
cuTexObjectGetResourceViewDesc:
    jmp *kt_real_cuTexObjectGetResourceViewDesc(%rip)
    .globl cuTexObjectGetTextureDesc
    .type cuTexObjectGetTextureDesc,@function
cuTexObjectGetTextureDesc:
    jmp *kt_real_cuTexObjectGetTextureDesc(%rip)
    .globl cuTexRefCreate
    .type cuTexRefCreate,@function
cuTexRefCreate:
    jmp *kt_real_cuTexRefCreate(%rip)
    .globl cuTexRefDestroy
    .type cuTexRefDestroy,@function
cuTexRefDestroy:
    jmp *kt_real_cuTexRefDestroy(%rip)
    .globl cuTexRefGetAddress
    .type cuTexRefGetAddress,@function
cuTexRefGetAddress:
    jmp *kt_real_cuTexRefGetAddress(%rip)
    .globl cuTexRefGetAddressMode
    .type cuTexRefGetAddressMode,@function
cuTexRefGetAddressMode:
    jmp *kt_real_cuTexRefGetAddressMode(%rip)
    .globl cuTexRefGetAddress_v2
    .type cuTexRefGetAddress_v2,@function
cuTexRefGetAddress_v2:
    jmp *kt_real_cuTexRefGetAddress_v2(%rip)
    .globl cuTexRefGetArray
    .type cuTexRefGetArray,@function
cuTexRefGetArray:
    jmp *kt_real_cuTexRefGetArray(%rip)
    .globl cuTexRefGetBorderColor
    .type cuTexRefGetBorderColor,@function
cuTexRefGetBorderColor:
    jmp *kt_real_cuTexRefGetBorderColor(%rip)
    .globl cuTexRefGetFilterMode
    .type cuTexRefGetFilterMode,@function
cuTexRefGetFilterMode:
    jmp *kt_real_cuTexRefGetFilterMode(%rip)
    .globl cuTexRefGetFlags
    .type cuTexRefGetFlags,@function
cuTexRefGetFlags:
    jmp *kt_real_cuTexRefGetFlags(%rip)
    .globl cuTexRefGetFormat
    .type cuTexRefGetFormat,@function
cuTexRefGetFormat:
    jmp *kt_real_cuTexRefGetFormat(%rip)
    .globl cuTexRefGetMaxAnisotropy
    .type cuTexRefGetMaxAnisotropy,@function
cuTexRefGetMaxAnisotropy:
    jmp *kt_real_cuTexRefGetMaxAnisotropy(%rip)
    .globl cuTexRefGetMipmapFilterMode
    .type cuTexRefGetMipmapFilterMode,@function
cuTexRefGetMipmapFilterMode:
    jmp *kt_real_cuTexRefGetMipmapFilterMode(%rip)
    .globl cuTexRefGetMipmapLevelBias
    .type cuTexRefGetMipmapLevelBias,@function
cuTexRefGetMipmapLevelBias:
    jmp *kt_real_cuTexRefGetMipmapLevelBias(%rip)
    .globl cuTexRefGetMipmapLevelClamp
    .type cuTexRefGetMipmapLevelClamp,@function
cuTexRefGetMipmapLevelClamp:
    jmp *kt_real_cuTexRefGetMipmapLevelClamp(%rip)
    .globl cuTexRefGetMipmappedArray
    .type cuTexRefGetMipmappedArray,@function
cuTexRefGetMipmappedArray:
    jmp *kt_real_cuTexRefGetMipmappedArray(%rip)
    .globl cuTexRefSetAddress
    .type cuTexRefSetAddress,@function
cuTexRefSetAddress:
    jmp *kt_real_cuTexRefSetAddress(%rip)
    .globl cuTexRefSetAddress2D
    .type cuTexRefSetAddress2D,@function
cuTexRefSetAddress2D:
    jmp *kt_real_cuTexRefSetAddress2D(%rip)
    .globl cuTexRefSetAddress2D_v2
    .type cuTexRefSetAddress2D_v2,@function
cuTexRefSetAddress2D_v2:
    jmp *kt_real_cuTexRefSetAddress2D_v2(%rip)
    .globl cuTexRefSetAddress2D_v3
    .type cuTexRefSetAddress2D_v3,@function
cuTexRefSetAddress2D_v3:
    jmp *kt_real_cuTexRefSetAddress2D_v3(%rip)
    .globl cuTexRefSetAddressMode
    .type cuTexRefSetAddressMode,@function
cuTexRefSetAddressMode:
    jmp *kt_real_cuTexRefSetAddressMode(%rip)
    .globl cuTexRefSetAddress_v2
    .type cuTexRefSetAddress_v2,@function
cuTexRefSetAddress_v2:
    jmp *kt_real_cuTexRefSetAddress_v2(%rip)
    .globl cuTexRefSetArray
    .type cuTexRefSetArray,@function
cuTexRefSetArray:
    jmp *kt_real_cuTexRefSetArray(%rip)
    .globl cuTexRefSetBorderColor
    .type cuTexRefSetBorderColor,@function
cuTexRefSetBorderColor:
    jmp *kt_real_cuTexRefSetBorderColor(%rip)
    .globl cuTexRefSetFilterMode
    .type cuTexRefSetFilterMode,@function
cuTexRefSetFilterMode:
    jmp *kt_real_cuTexRefSetFilterMode(%rip)
    .globl cuTexRefSetFlags
    .type cuTexRefSetFlags,@function
cuTexRefSetFlags:
    jmp *kt_real_cuTexRefSetFlags(%rip)
    .globl cuTexRefSetFormat
    .type cuTexRefSetFormat,@function
cuTexRefSetFormat:
    jmp *kt_real_cuTexRefSetFormat(%rip)
    .globl cuTexRefSetMaxAnisotropy
    .type cuTexRefSetMaxAnisotropy,@function
cuTexRefSetMaxAnisotropy:
    jmp *kt_real_cuTexRefSetMaxAnisotropy(%rip)
    .globl cuTexRefSetMipmapFilterMode
    .type cuTexRefSetMipmapFilterMode,@function
cuTexRefSetMipmapFilterMode:
    jmp *kt_real_cuTexRefSetMipmapFilterMode(%rip)
    .globl cuTexRefSetMipmapLevelBias
    .type cuTexRefSetMipmapLevelBias,@function
cuTexRefSetMipmapLevelBias:
    jmp *kt_real_cuTexRefSetMipmapLevelBias(%rip)
    .globl cuTexRefSetMipmapLevelClamp
    .type cuTexRefSetMipmapLevelClamp,@function
cuTexRefSetMipmapLevelClamp:
    jmp *kt_real_cuTexRefSetMipmapLevelClamp(%rip)
    .globl cuTexRefSetMipmappedArray
    .type cuTexRefSetMipmappedArray,@function
cuTexRefSetMipmappedArray:
    jmp *kt_real_cuTexRefSetMipmappedArray(%rip)
    .globl cuThreadExchangeStreamCaptureMode
    .type cuThreadExchangeStreamCaptureMode,@function
cuThreadExchangeStreamCaptureMode:
    jmp *kt_real_cuThreadExchangeStreamCaptureMode(%rip)
    .globl cuUserObjectCreate
    .type cuUserObjectCreate,@function
cuUserObjectCreate:
    jmp *kt_real_cuUserObjectCreate(%rip)
    .globl cuUserObjectRelease
    .type cuUserObjectRelease,@function
cuUserObjectRelease:
    jmp *kt_real_cuUserObjectRelease(%rip)
    .globl cuUserObjectRetain
    .type cuUserObjectRetain,@function
cuUserObjectRetain:
    jmp *kt_real_cuUserObjectRetain(%rip)
    .globl cuVDPAUCtxCreate
    .type cuVDPAUCtxCreate,@function
cuVDPAUCtxCreate:
    jmp *kt_real_cuVDPAUCtxCreate(%rip)
    .globl cuVDPAUCtxCreate_v2
    .type cuVDPAUCtxCreate_v2,@function
cuVDPAUCtxCreate_v2:
    jmp *kt_real_cuVDPAUCtxCreate_v2(%rip)
    .globl cuVDPAUGetDevice
    .type cuVDPAUGetDevice,@function
cuVDPAUGetDevice:
    jmp *kt_real_cuVDPAUGetDevice(%rip)
    .globl cuWaitExternalSemaphoresAsync
    .type cuWaitExternalSemaphoresAsync,@function
cuWaitExternalSemaphoresAsync:
    jmp *kt_real_cuWaitExternalSemaphoresAsync(%rip)
    .globl cuWaitExternalSemaphoresAsync_ptsz
    .type cuWaitExternalSemaphoresAsync_ptsz,@function
cuWaitExternalSemaphoresAsync_ptsz:
    jmp *kt_real_cuWaitExternalSemaphoresAsync_ptsz(%rip)
    .section .note.GNU-stack,"",@progbits
