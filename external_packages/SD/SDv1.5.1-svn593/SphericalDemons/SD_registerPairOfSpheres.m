function parms = SD_registerPairOfSpheres(fixed_sphere, moving_sphere, SD_atlas, reg_parms)

% SD_registerPairOfSpheres(fixed_sphere, moving_sphere, SD_atlas, parms)
%
% Diffeomorphic registration of 2 spheres.
% Calls SD_registerAtlas2SphereMultiRes by creating a fake atlas using moving_sphere. 
%
% Just need to pass in a basic SD_atlas shell. 

num_meshes = length(SD_atlas.parms.uniform_meshes);

for i = num_meshes:-1:1
    disp(['Computing basic atlas ' num2str(i)]);
    basic_atlas = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, SD_atlas.parms.uniform_meshes{i});
    
    reg_parms.meshes{i} = basic_atlas;
    
    basic_atlas.vertices = MARS_yrotate(MARS_zrotate(basic_atlas.vertices, SD_atlas.zrotate), SD_atlas.yrotate);

    if(SD_atlas.multidim)
        if(i == num_meshes)
           basic_atlas.mean = MARS_linearInterpolate(basic_atlas.vertices, moving_sphere, moving_sphere.data);
        else
           temp             = MARS_AverageData(SD_atlas.basic_atlas{i+1}, SD_atlas.basic_atlas{i+1}.mean, 0, 1);
           basic_atlas.mean = temp(:, 1:size(basic_atlas.vertices, 2));
        end
    else   
        if(i == num_meshes)
            temp = MARS_linearInterpolate(basic_atlas.vertices, moving_sphere, moving_sphere.data);
        else
            temp = MARS_AverageData(SD_atlas.basic_atlas{i+1}, temp, 0, 1);
        end
    
	temp = temp(:, 1:size(basic_atlas.vertices, 2));
        basic_atlas.mean = temp(i, :);    
    end
    
    basic_atlas.variance = zeros(size(basic_atlas.mean), 'single') + 1;
    SD_atlas.basic_atlas{i} = basic_atlas;
end
parms = SD_registerAtlas2SphereMultiRes(fixed_sphere, SD_atlas, reg_parms);


