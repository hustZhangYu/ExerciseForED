% Regress the actual periodic Ising exercise against independent spin energies.
testRoot = fileparts(fileparts(mfilename('fullpath')));
exerciseSource = fileread(fullfile(testRoot, 'code1.txt'));
exerciseBlock = regexp(exerciseSource, ...
    '(?s)%Exercise2.*?(?=%%\s*%Exercise3)', 'match', 'once');
assert(~isempty(exerciseBlock), 'Ising exercise section not found.');
eval(exerciseBlock);
% The exercise's clear-all also clears the test setup; its L, J, basis and Ham remain.
spins = 2*basis-1;
expected = J/4 * sum(spins .* circshift(spins, [0, 1]), 2);
assert(norm(diag(Ham)-expected) < 1e-12, ...
    'Periodic Ising diagonal must contain every nearest-neighbour bond.');
disp('Periodic Ising exercise check passed.');

