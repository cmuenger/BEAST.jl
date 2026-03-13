
# Configure and launch a GPU kernel, given the kernel function and its arguments.
function kernel_config!(gpu_kernel,args...)
    kernel = @cuda launch=false gpu_kernel(args...)
    config = launch_configuration(kernel.fun)
    
    threads = config.threads
    blocks  = config.blocks

    CUDA.@sync begin
         kernel(args...; threads, blocks)
    end

    return 
end