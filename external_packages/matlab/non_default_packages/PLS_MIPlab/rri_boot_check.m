function [min_subj_per_group, is_boot_samples, boot_samples, new_num_boot] ...
	= rri_boot_check(num_subj_lst, num_cond, num_boot, incl_seq)

   if ~exist('incl_seq','var')
      incl_seq = 0;
   end

   total_subj = sum(num_subj_lst);
   num_group = length(num_subj_lst);
   total_rows = num_cond * total_subj;

   min_subj_per_group = [];
   is_boot_samples = [];
   boot_samples = [];

   new_num_boot = num_boot;
   done = 0;

   %  num_subj in one of group is less than 3
   %
   if (min(num_subj_lst) < 3)
      error('Bootstrap analysis requires that each group must have at least 3 subjects');
   end;

   percentage = 50;

   min_subj_per_group = ceil(min(num_subj_lst)*percentage/100);
   max_subj_per_group = 8;

   %  if num of subj is 8 or less, we can get from boot matrix
   %
   is_boot_samples = zeros(1,num_group);
   boot_samples = cell(1,num_group);

   vernum = get_matlab_version;
   if vernum < 7002
      rand('state',sum(100*clock));
   elseif vernum < 7013
      rand('twister',sum(100*clock));
   else
      rng_default; rng_shuffle;
   end

   if (sum(num_subj_lst <= max_subj_per_group) == num_group)

      for g = 1:num_group

         num_subj = num_subj_lst(g);
%         diff_subj = min_subj_per_group;

         boot_sample2 = [];

         %  remove sequential order
         %
         if incl_seq
            max_diff_subj = num_subj;
         else
            max_diff_subj = num_subj - 1;
         end

         for diff_subj = min_subj_per_group:max_diff_subj
            boot_sample1 = rri_boot_samples(num_subj, diff_subj);
            boot_sample2 = [boot_sample2; boot_sample1];
         end

         num_boot_samples = size(boot_sample2, 1);

         if num_boot_samples < new_num_boot
            disp(['WARNING:' num2str(num_subj) ' subjects can only have ' num2str(num_boot_samples) ' different bootstrap samples. Set num_boot to this theoretical maximum.']);
            new_num_boot = num_boot_samples;
         end

         if ~isempty(boot_sample2)
            boot_sample2 = boot_sample2(randperm(num_boot_samples),:);
            boot_samples{g} = boot_sample2;
            is_boot_samples(g) = 1;
         end

      end	% for
   end		% if

   return;

