clc;
close all;
    
crc2 = [1 1 1];
crc4 = [1 1 0 1];
crc6 = [1 1 0 1 1 1 1];
crc8 = [1 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 1];
crc1 = [1 1];
    Pe_decode1;
    P_ed1;
    Pe_decode2;
    P_ed2 ;
    Pe_decode3;
    P_ed3 ;
    Pe_decode4;
    P_ed4 ;
    Pe_decode5;
    P_ed5;
M = [2, 4, 6, 8, 1]; %����� ������������������ 

for ih = 1 : 5
    crc;
    m = 0;
if ih == 1
    m = M(ih);
    crc = crc2;
end
if ih == 2
    m = M(ih);
    crc = crc4;
end
if ih == 3
    m = M(ih);
    crc = crc6;
end
if ih == 4
    m = M(ih);
    crc = crc8;
end
if ih == 5
    m = M(ih);
    crc = crc1;
end

r = length(crc)-1; %2^r ����� ��������� ����
n = 1000; %���������� �������
%CNR ��� ������� ����� ������ �������� �� -inf �� 1 ��� ���
snr_arr = -30:5:10;
snr_theor_arr = -30:10;
%������������� �������� ����������� ������ �� ��� �������� �� ������� ��
%������ q(sqrt(2E/N0)) %��� ��������� SNR
Pb_theor =  qfunc(sqrt(2*(10.^(snr_theor_arr./10))));

%������������� �������� ����������� ������ ������������� (1/2)^r
Pe_decode_theor = ones(1,length(snr_theor_arr)).*(1/2^r);

%�� � ������ ������
%��������� ���� ������� ����
codewords = zeros(2^m, m+r);
%��� BPSK
modulate_codewords = zeros(2^m, m+r);
d_arr = zeros((2^m), 1);

for word = 0:bi2de(ones(1, m))%���������� ��������� ������
    %��������� � ����������
    %������� ��������� � �������� ����
    [~, c] = gfdeconv(de2bi(bitshift(word, r)), crc);
    %����������� ����� � ������ ����
    a = de2bi(bitxor((bitshift(word, r)), bi2de(c)), m+r);
    %�������� ������ �����
    codewords(word+1, :) = a;
    %�������� ���������� ��������� ��������� ������� �������
    d_arr(word+1) = nnz(a);
    %������� �� 1 � -1 BPSK
    a_m = de2bi(a).*(-2)+1; 
    modulate_codewords(word+1, :) = a_m;
end
%��������� ������ ����������� ������ �������������
d_min = min(d_arr(2:end));
P_ed = zeros(1, length(Pb_theor));
index = 1;
for p = Pb_theor
    tmp = 0;
    for i = d_min:(m+r)
        Ai = sum(d_arr == i);
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i));
    end
    P_ed(index) = tmp;
    index = index + 1;
end

index = 1;
%���� ��� ��������� SNR
for SNR = snr_arr
    N_b = 0;
    Ne_decode = 0;
    %��� ������� ������
    for i = 1:n
        %��������� ���������
        rnd_ind = randi(2^m,1);
        massage = modulate_codewords(rnd_ind, :);
        %��������� � ����� � �������� ����
        SNRi = 10^(SNR/10);
        sigma = sqrt(1/(2*SNRi));
        b_m = massage+sigma*randn(1, length(a_m));
        %������ �� ���� ������������
        b = b_m < 0;
        %� ��� �������
        a = codewords(rnd_ind, :); %����� ��������� ���������
        e = xor(a,b);%�������� ���������� �������� ������ ������
        N_b = N_b + nnz(e); %����� ������ �� ��� 
        %����� ������ �������������
        [~, s] = gfdeconv(double(b), crc);%����� ���� ������� 0 �� ��� ��
         if (bi2de(s) == 0) & (nnz(e) ~= 0)
            Ne_decode = Ne_decode + 1;
         end
    end
    
    Pb(index) = N_b/(n*(m+r));
    Pe_decode(index) = Ne_decode/n;
    index = index + 1;
end
  
if ih == 1
   Pe_decode1 =  Pe_decode;
   P_ed1 = P_ed;
end
if ih == 2
    Pe_decode2 =  Pe_decode;
    P_ed2 = P_ed;
end
if ih == 3
   Pe_decode3 =  Pe_decode;
   P_ed3 = P_ed;
end
if ih == 4
    Pe_decode4 =  Pe_decode;
    P_ed4 = P_ed;
end
if ih == 5
    Pe_decode5 =  Pe_decode;
    P_ed5 = P_ed;
end


end
%������ ������� 

figure(2);

semilogy(snr_arr, Pe_decode1, "bo", snr_theor_arr, P_ed1, snr_arr, Pe_decode2, "bo", snr_theor_arr, P_ed2, snr_arr, Pe_decode3, "bo", snr_theor_arr, P_ed3, snr_arr, Pe_decode4, "bo", snr_theor_arr, P_ed4, snr_arr, Pe_decode5, "bo", snr_theor_arr, P_ed5);




