function hf = CombFigs(varargin)
% CombFigs��varargin������ͬfiguresͼƬ�ϲ���һ��figure��
% ���ø�ʽ 
% ����ʽ�� 
%          CombFigs()  �����κ��������������ǰĿ¼������fig�ļ��ϲ���һ��fig�ļ�
%                      Ĭ�Ϻϲ����fig�ļ���Ϊ  �ϲ�Figure�ļ�.fig
% ָ��Ŀ��fig�ļ�����
%          CombFigs(figname)  fignameΪ�ϲ����fig�ļ���������ǰĿ¼������fig�ļ��ϲ�Ϊfigname�ļ� 
% ָ��Ŀ��fig�ļ�����Ҫ�ϲ����ļ�����
%          CombFigs(figname��figfiles)  fignameΪ�ϲ����fig�ļ���,
%                                       figfilesΪ���ϲ���fig�ļ���cell�����ʽ����figfiles�ļ��ϲ�Ϊfigname�ļ� 
% By ZFS@wust 20210905
% ��ȡ����Matlab/Simulinkԭ�����Ϻͳ������ע΢�Ź��ںţ�Matlab Fans

if  isempty(varargin)                   % ������
    files = dir( '*.fig' );
    filenames = { files(:).name };
    filenames0 = filenames{1};        % �ϲ�Ŀ���ļ�
    figname = '�ϲ�Figure�ļ�.fig';    % Ĭ�ϱ����ļ���
elseif length(varargin) == 1          % ����һ��Ԫ��ʱ��Ϊ����Ŀ���ļ���
    files = dir('*.fig');
    filenames = {files(:).name};
    filenames0 = filenames{1};     % �ϲ�Ŀ���ļ�
    figname = varargin{1};         % �����ļ���
else                               % ������Ԫ��ʱ��ÿ��Ԫ��Ϊfig�ļ���
    figname = varargin{1};       % ��1�������������ļ���
    filenames  = varargin{2};      % ��2�����������ϲ���fig�ļ��б��ַ���cell����
    filenames0 = filenames{1};     % �ϲ�Ŀ���ļ�
end

if isempty(figname)                  % ��filenames2����Ϊ��ռλ��ʱ
    figname = '�ϲ�Figure�ļ�.fig';   % Ĭ�ϱ����ļ���
end

if isempty(filenames)
    error('û�пɹ��ϲ���figures�ļ�')
else
    hf = open(filenames0);
    ax = findall(hf,'type','axes');      % axes���
    hg = findall(hf,'type','legend');   % legend���
    if isempty(hg)
        for ii=1:length(allchild(ax))
            hS(ii)={[filenames0  '����' num2str(ii)]};
        end
    else
        hS=hg.String;
    end
end

for ii = 1:length( filenames )
    if strcmp(filenames{ii},filenames0)
        continue;
    end
    h2= openfig(filenames{ii}, 'invisible');
    ax2 = findall(h2,'type','axes');
    hL = allchild(ax2);
    copyobj(hL,ax);   % ��������
    hg2 = findall(h2,'type','legend');
    if isempty(hg2)
        for kk = 1:length(hL)
            jj = length(hS);
            hS(jj+1) = { [filenames{ii} '����' num2str(kk)]};
        end
    else
        jj = length(hS);
        n = length(hg2.String);
        hS(jj+1:jj+n) = hg2. String;
    end
end
%legend(ax,hS);
legend off
%savefig(hf,figname)          % ����figure�ļ�

end


