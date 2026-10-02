function eventInds = findEvents_dyn(regime)
% Start and end indices of the fully dynamic phases of a
% quasidynamic_and_dynamic run.
%
% regime is the saved /time/regime1D (or regime2D) series: 0 for
% quasi-dynamic steps, 1 for fully dynamic steps. eventInds has one row per
% event, [first dynamic index, last dynamic index]; an event still running
% at the end of the output ends at the last index.

regime = regime(:) > 0.5;
d = diff([false; regime; false]);
startInds = find(d == 1);
endInds = find(d == -1) - 1;
eventInds = [startInds endInds];
end
