clc;
close all;
%crc = [1 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 1];  %crc16 m=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
crc = [1 0 1 0 1 0 1 1 1]; %crc8 m=4 x8+x7+x6+x4+x2+1
m = 4; %длина последовательности 
r = length(crc)-1; %2^r колво элементов поля
n = 100000; %количество пакетов
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
%строим графики 
figure(1);
semilogy(snr_theor_arr, Pb_theor, 'c', snr_arr, Pb, 'ko');
grid on
figure(2);
semilogy(snr_theor_arr, Pe_decode_theor, 'm', snr_arr, Pe_decode, 'ko', snr_theor_arr, P_ed, 'c--');
grid on





