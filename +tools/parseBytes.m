function [outputStr] = parseBytes(str)
%UNTITLED3 Function to parse unicode bytestring and generate display character

split_str  = split(str, "\");
split_str  = split_str(2:end);
split_str = cellfun(@(x) x(2:end), split_str, 'UniformOutput',false);
a = cell2mat(split_str);
b = a';
c = b(:);
d = c';
n = length(split_str);
bin_rep = hexToBinaryVector(d, n*8);

switch n
    case 1
        final_hex = binaryVectorToHex(bin_rep(2:8));
    case 2
        final_hex = binaryVectorToHex([bin_rep(4:8) bin_rep(11:16)]);
    case 3
        final_hex = binaryVectorToHex([bin_rep(5:8) bin_rep(11:16) bin_rep(19:24)]);
    case 4
        final_hex = binaryVectorToHex([bin_rep(6:8) bin_rep(11:16) bin_rep(19:24) bin_rep(27:32)]);
        
end
outputStr = char(hex2dec(final_hex));
end