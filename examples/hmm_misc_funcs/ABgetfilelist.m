function filelist = ABgetfilelist(dirname,filestr)
% filelist = ABgetfilelist(dirname,filestr)
    filelist = ls(fullfile(dirname,filestr));
    filelist = textscan(filelist,'%[^\n]');
    filelist = filelist{1};
    filelist = sort_nat(filelist);
end