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
M = [2, 4, 6, 8, 1]; %длина последовательности 

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

r = length(crc)-1; %2^r колво элементов поля
n = 1000; %количество пакетов
%CNR при которых будем делать проверку от -inf до 1 для Дец
snr_arr = -30:5:10;
snr_theor_arr = -30:10;
%теоретические значения вероятности ошибки на бит вычислим по формуле из
%методы q(sqrt(2E/N0)) %для различных SNR
Pb_theor =  qfunc(sqrt(2*(10.^(snr_theor_arr./10))));

%теоретическое значение вероятности ошибки декодирования (1/2)^r
Pe_decode_theor = ones(1,length(snr_theor_arr)).*(1/2^r);

%ну а теперь основа
%полученеи всех кодовых слов
codewords = zeros(2^m, m+r);
%для BPSK
modulate_codewords = zeros(2^m, m+r);
d_arr = zeros((2^m), 1);

for word = 0:bi2de(ones(1, m))%модулирует случайный сигнал
    %проблемки с полиномами
    %деление полиномов в конечном поле
    [~, c] = gfdeconv(de2bi(bitshift(word, r)), crc);
    %преобразуем числа в екторы цифр
    a = de2bi(bitxor((bitshift(word, r)), bi2de(c)), m+r);
    %получили кодове слово
    codewords(word+1, :) = a;
    %получили количество ненулевых элементов энергия сигнала
    d_arr(word+1) = nnz(a);
    %сделали из 1 и -1 BPSK
    a_m = de2bi(a).*(-2)+1; 
    modulate_codewords(word+1, :) = a_m;
end
%посчитаем точную вероятность ошибки декодирования
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
%цикл для различных SNR
for SNR = snr_arr
    N_b = 0;
    Ne_decode = 0;
    %для каждого пакета
    for i = 1:n
        %появилось сообщение
        rnd_ind = randi(2^m,1);
        massage = modulate_codewords(rnd_ind, :);
        %отправили в канал и наложили АБГШ
        SNRi = 10^(SNR/10);
        sigma = sqrt(1/(2*SNRi));
        b_m = massage+sigma*randn(1, length(a_m));
        %пришли на вход демодулятора
        b = b_m < 0;
        %а тут декодер
        a = codewords(rnd_ind, :); %берем известное сообщение
        e = xor(a,b);%побитово складываем получаем вектор ошибок
        N_b = N_b + nnz(e); %колво ошибок на бит 
        %колво ошибок декодирования
        [~, s] = gfdeconv(double(b), crc);%делим если остаток 0 то все ок
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
%строим графики 

figure(2);

semilogy(snr_arr, Pe_decode1, "bo", snr_theor_arr, P_ed1, snr_arr, Pe_decode2, "bo", snr_theor_arr, P_ed2, snr_arr, Pe_decode3, "bo", snr_theor_arr, P_ed3, snr_arr, Pe_decode4, "bo", snr_theor_arr, P_ed4, snr_arr, Pe_decode5, "bo", snr_theor_arr, P_ed5);




