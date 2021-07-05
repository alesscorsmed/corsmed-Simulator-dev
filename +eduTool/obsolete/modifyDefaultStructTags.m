function struct_tags = modifyDefaultStructTags(struct_tags)

if ~isfield(struct_tags,'CUDA')
    struct_tags.CUDA = 8;
end

if ~isfield(struct_tags,'parfeval')
    struct_tags.parfeval = 0;
end

if ~isfield(struct_tags,'PYTHON')
    struct_tags.PYTHON = '3.7';
end

if ~isfield(struct_tags,'Development')
    struct_tags.Development = '0';
end

if ~isfield(struct_tags,'activateNumerical')
    struct_tags.activateNumerical = '0';
end