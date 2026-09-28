sig_levels = NaN(length(vbl_old.dirList), 1);
image_names = strings(length(vbl.dirList), 1);

for d = 1:length(vbl_old.dirList)
    sig_levels(d) = compare_real_vs_rand2(vbl_old, d);
    image_names(d) = string(vbl_old.dirList(d).name);
end

save('all_sig_levelsOLD.mat', 'sig_levels', 'image_names');