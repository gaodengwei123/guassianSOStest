function display(E)

fprintf('\n');
disp([inputname(1) ' =']);
[m, ~] = size(E.X);

fprintf('%d obstacle in Map\n', m);

fprintf('\n');
fprintf('range of Map:\n'); disp(E.range);
fprintf('vertex of Map:\n'); disp(E.vertex);

if isempty(E)
    fprintf('Empty Map.\n\n');
    return;
end

return;
