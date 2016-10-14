function r = CBIG_find_maxLikeRun(kDir)

dirList = dir([kDir 'r*']);

maxLogLike = -inf;

figure;
hold on;

for idx = 1:numel(dirList)
    runName = dirList(idx).name;
    try
        logLikes = load([kDir runName '/likelihood.dat']);
        logLike = logLikes(end, 1);
        runID = str2num(runName(2:end));
        scatter(runID, logLike, 'b.');
        if logLike > maxLogLike
            maxLogLike = logLike;
            maxRunID = runID;
        end
    catch
        warning([kDir runName ': no likelihood.dat']);
    end
end

scatter(maxRunID, maxLogLike, 'ro');
hold off;
ylabel('Log-likelihood');
xlabel('Run');

r = num2str(maxRunID);
