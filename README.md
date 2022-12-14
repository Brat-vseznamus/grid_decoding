# Построение решетки по порождающей матрице

Чтобы получить решетку нужно указать в `matrix.txt` порождающую матрицу, далее она сведётся с МСФ, и граф решетки будет построен в `graph.txt` в Graphviz аннотации. Чтобы узнать какие на каждом слое используются активные символы и в каком они порядке запустите:

```py
index = #нужный столбец

print(active_rows(index, rngs))
```

Чтобы восстановить Y, укажите его в файле `y.txt` в таком же формате, что и матрицу G, в конце будет выведен оригинал и исправленная версия.

## Пример вывода

Данные взяты из прикрепленных файлов

```
МФС:
[[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0],
 [0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1],
 [1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
 [0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0],
 [0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0]]

log(layer[0]) = 0
log(layer[1]) = 1
log(layer[2]) = 2
log(layer[3]) = 3
log(layer[4]) = 4
log(layer[5]) = 5
log(layer[6]) = 5
log(layer[7]) = 6
log(layer[8]) = 6
log(layer[9]) = 6
log(layer[10]) = 6
log(layer[11]) = 5
log(layer[12]) = 4
log(layer[13]) = 3
log(layer[14]) = 2
log(layer[15]) = 1
log(layer[16]) = 0

original: [0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1]
fixed:    [0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1]
```