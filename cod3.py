#######################################################
##     TRABALHO FINAL DE COMPUTAÇÃO GRÁFICA 2025     ##
## ALUNOS: LUIS EDUARDO FURTADO CONDE (202111140024) ##   
## RICARDO MARIANO NUNES CORTINHAS (202011140025)    ##
#######################################################

import tkinter as tk
import numpy as np
import sys
from tkinter import simpledialog, messagebox

class Pratica_CG:
    def __init__(self, root):
        self.root = root
        self.root.title("Trabalho de Computação Gráfica")
        
        # Configurações do Canvas e Escala 
        self.canvas_width = 999
        self.canvas_height = 500
        self.canvas = tk.Canvas(root, width=self.canvas_width, height=self.canvas_height, bg="#f0f0f0")
        self.canvas.pack(pady=10)
        
        # Fator de escala para converter coordenadas do plano para pixels
        self.escala = 20
        self.center_x = self.canvas_width / 2
        self.center_y = self.canvas_height / 2

        # matriz de pixels que armazena o estado da imagem
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)
        
        self.desenhar_plano_cartesiano()
        
        # Variável para armazenar os vértices  ativos
        self.vertices_ativos = None
        self.forma_ativa_tag = "forma"
        

        # Frames para a organização dos botões 
        self.frame_botoes_forma = tk.Frame(root)
        self.frame_botoes_forma.pack(pady=5)
        
        tk.Label(self.frame_botoes_forma, text="Desenhar Forma:").pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Reta (Bresenham)", command=self.solicitar_reta).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Círculo", command=self.solicitar_circulo).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Curva de Bézier", command=self.solicitar_bezier).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Polilinha", command=self.solicitar_polilinha).pack(side="left", padx=5)

        self.frame_botoes_transf = tk.Frame(root)
        self.frame_botoes_transf.pack(pady=5)
        
        tk.Label(self.frame_botoes_transf, text="Transformações:").pack(side="left", padx=5)
        tk.Button(self.frame_botoes_transf, text="Translação", command=self.solicitar_transladar).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_transf, text="Rotação", command=self.solicitar_rotacionar).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_transf, text="Escala", command=self.solicitar_escalar).pack(side="left", padx=5)

        self.frame_processamento = tk.Frame(root)
        self.frame_processamento.pack(pady=5)
        
        tk.Label(self.frame_processamento, text="Processamento:").pack(side="left", padx=5)
        tk.Button(self.frame_processamento, text="Preenchimento Recursivo", command=self.solicitar_preenchimento_recursivo).pack(side="left", padx=5)
        tk.Button(self.frame_processamento, text="Preenchimento Varredura", command=self.solicitar_preenchimento_varredura).pack(side="left", padx=5)
        tk.Button(self.frame_processamento, text="Recorte", command=self.solicitar_recorte).pack(side="left", padx=5)
        
        self.btn_limpar = tk.Button(root, text="Limpar", command=self.limpar_canvas)
        self.btn_limpar.pack(pady=10)

    #  Funções do Plano Cartesiano 
    def desenhar_plano_cartesiano(self):
        """Desenha a grade e os eixos X e Y no canvas."""
        self.canvas.delete("grid")

        # Desenhar linhas verticais
        for i in range(0, self.canvas_width, self.escala):
            self.canvas.create_line(i, 0, i, self.canvas_height, fill="#ddd", tags="grid")
        # Desenhar linhas horizontais
        for i in range(0, self.canvas_height, self.escala):
            self.canvas.create_line(0, i, self.canvas_width, i, fill="#ddd", tags="grid")

        # Desenhar eixos principais
        self.canvas.create_line(0, self.center_y, self.canvas_width, self.center_y, fill="#333", width=2, tags="grid")
        self.canvas.create_line(self.center_x, 0, self.center_x, self.canvas_height, fill="#333", width=2, tags="grid")
        
    def converter_para_canvas(self, x, y):
        """Converte coordenadas do plano para pixels do canvas."""
        return self.center_x + x * self.escala, self.center_y - y * self.escala

    def desenhar_grid_no_canvas(self):
        """Desenha o conteúdo da matriz de pixels no canvas."""
        self.canvas.delete(self.forma_ativa_tag)
        cores = {0: "", 1: "blue", 2: "red"}
        for x in range(self.canvas_width):
            for y in range(self.canvas_height):
                cor_valor = self.pixel_grid[x, y]
                if cor_valor in cores and cores[cor_valor] != "":
                    self.canvas.create_rectangle(x, y, x + 1, y + 1, fill=cores[cor_valor], outline=cores[cor_valor], tags=self.forma_ativa_tag)
    
    #  Funções para Desenho de Formas 
    def limpar_formas(self):
        """Limpa apenas as formas desenhadas, mantendo a grade."""
        self.canvas.delete(self.forma_ativa_tag)
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)

    def desenhar_reta_bresenham(self, x0, y0, x1, y1):
        """Implementa o algoritmo de Bresenham para desenhar uma linha."""
                       
        x0_p, y0_p = self.converter_para_canvas(x0, y0)
        x1_p, y1_p = self.converter_para_canvas(x1, y1)
        dx = abs(x1_p - x0_p)
        dy = abs(y1_p - y0_p)
        sx = 1 if x0_p < x1_p else -1
        sy = 1 if y0_p < y1_p else -1
        err = dx - dy
        
        x_p, y_p = int(x0_p), int(y0_p)
        while True:
            if 0 <= x_p < self.canvas_width and 0 <= y_p < self.canvas_height:
                self.pixel_grid[x_p, y_p] = 1
            if x_p == int(x1_p) and y_p == int(y1_p): break
            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                x_p += sx
            if e2 < dx:
                err += dx
                y_p += sy
        
        self.vertices_ativos = np.array([[x0, y0], [x1, y1]])
        self.desenhar_grid_no_canvas()
 
    def desenhar_poligono(self, vertices):
        """Desenha um polígono a partir de uma lista de vértices [x, y]."""
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)
        self.limpar_formas()
        # Converte as coordenadas do plano cartesiano para pixels do canvas
        pontos = []
        for x, y in vertices:
            pontos.append(self.center_x + x * self.escala)
            pontos.append(self.center_y - y * self.escala)
        self.canvas.create_polygon(pontos, outline="blue", fill="", width=3, tags=self.forma_ativa_tag)
        self.vertices_ativos = vertices

    def desenhar_poligono2(self, vertices):
        """Desenha um polígono na matriz de pixels."""
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)
        self.limpar_formas()
        for i in range(len(vertices)):
            x0, y0 = vertices[i]
            x1, y1 = vertices[(i + 1) % len(vertices)]
            self.desenhar_reta_bresenham_interno(x0, y0, x1, y1)
        
        self.desenhar_grid_no_canvas()
        self.vertices_ativos = vertices

    def desenhar_reta_bresenham_interno(self, x0, y0, x1, y1):
        x0_p, y0_p = self.converter_para_canvas(x0, y0)
        x1_p, y1_p = self.converter_para_canvas(x1, y1)
        dx = abs(x1_p - x0_p)
        dy = abs(y1_p - y0_p)
        sx = 1 if x0_p < x1_p else -1
        sy = 1 if y0_p < y1_p else -1
        err = dx - dy
        
        x_p, y_p = int(x0_p), int(y0_p)
        while True:
            if 0 <= x_p < self.canvas_width and 0 <= y_p < self.canvas_height:
                self.pixel_grid[x_p, y_p] = 1
            if x_p == int(x1_p) and y_p == int(y1_p): break
            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                x_p += sx
            if e2 < dx:
                err += dx
                y_p += sy
  
    def desenhar_circulo(self, cx, cy, raio):
        """Desenha um círculo a partir de centro e raio."""
        self.limpar_formas()
        # Converte as coordenadas para o canvas
        x1 = self.center_x + (cx - raio) * self.escala
        y1 = self.center_y - (cy + raio) * self.escala
        x2 = self.center_x + (cx + raio) * self.escala
        y2 = self.center_y - (cy - raio) * self.escala
        x_p, y_p = self.converter_para_canvas(cx, cy)
        raio_p = raio * self.escala
        self.vertices_ativos = np.array([[cx, cy]])

        x = int(raio_p)
        y = 0
        p = 1 - raio_p

        while y <= x:
            self.desenhar_simetria_circulo(x_p, y_p, x, y, 1) # 1 para borda
            y += 1
            if p <= 0:
                p = p + 2 * y + 1
            else:
                x -= 1
                p = p + 2 * y - 2 * x + 1

        self.desenhar_grid_no_canvas()

    def desenhar_simetria_circulo(self, cx, cy, x, y, cor):
        pixels = [(x, y), (-x, y), (x, -y), (-x, -y), (y, x), (-y, x), (y, -x), (-y, -x)]
        for dx, dy in pixels:
            px, py = int(cx + dx), int(cy + dy)
            if 0 <= px < self.canvas_width and 0 <= py < self.canvas_height:
                self.pixel_grid[px, py] = cor

    def desenhar_linha_pixel(self, x0, y0, x1, y1):
        """Desenha uma linha na matriz de pixels usando Bresenham."""
        x0_p, y0_p = self.converter_para_canvas(x0, y0)
        x1_p, y1_p = self.converter_para_canvas(x1, y1)
        
        x0_p, y0_p = int(x0_p), int(y0_p)
        x1_p, y1_p = int(x1_p), int(y1_p)
        
        dx = abs(x1_p - x0_p)
        dy = abs(y1_p - y0_p)
        sx = 1 if x0_p < x1_p else -1
        sy = 1 if y0_p < y1_p else -1
        err = dx - dy
        
        x_p, y_p = x0_p, y0_p
        while True:
            if 0 <= x_p < self.canvas_width and 0 <= y_p < self.canvas_height:
                self.pixel_grid[x_p, y_p] = 1
            if x_p == x1_p and y_p == y1_p: break
            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                x_p += sx
            if e2 < dx:
                err += dx
                y_p += sy

    def desenhar_bezier(self, pontos_controle):
        self.limpar_formas()
                
        # Número de segmentos menor para evitar travamentos
        num_segmentos = 50
        
        # Desenhar os pontos de controle como pequenos círculos para referência
        for px, py in pontos_controle:
            x_p, y_p = self.converter_para_canvas(px, py)
            x_p, y_p = int(x_p), int(y_p)
            for dx in range(-2, 3):
                for dy in range(-2, 3):
                    nx, ny = x_p + dx, y_p + dy
                    if 0 <= nx < self.canvas_width and 0 <= ny < self.canvas_height:
                        if dx*dx + dy*dy <= 4:  # Aproximação de um círculo pequeno
                            self.pixel_grid[nx, ny] = 2  # Em vermelho
        
        pontos_curva = []
        
        if len(pontos_controle) == 2:
            # Caso de linha reta
            pontos_curva = [pontos_controle[0], pontos_controle[1]]
        elif len(pontos_controle) == 3:
            # Curva quadrática: B(t) = (1-t)²P₀ + 2(1-t)tP₁ + t²P₂
            for i in range(num_segmentos + 1):
                t = i / num_segmentos
                x = (1-t)**2 * pontos_controle[0][0] + 2*(1-t)*t * pontos_controle[1][0] + t**2 * pontos_controle[2][0]
                y = (1-t)**2 * pontos_controle[0][1] + 2*(1-t)*t * pontos_controle[1][1] + t**2 * pontos_controle[2][1]
                pontos_curva.append((x, y))
        elif len(pontos_controle) == 4:
            # Curva cúbica: B(t) = (1-t)³P₀ + 3(1-t)²tP₁ + 3(1-t)t²P₂ + t³P₃
            for i in range(num_segmentos + 1):
                t = i / num_segmentos
                x = (1-t)**3 * pontos_controle[0][0] + 3*(1-t)**2*t * pontos_controle[1][0] + 3*(1-t)*t**2 * pontos_controle[2][0] + t**3 * pontos_controle[3][0]
                y = (1-t)**3 * pontos_controle[0][1] + 3*(1-t)**2*t * pontos_controle[1][1] + 3*(1-t)*t**2 * pontos_controle[2][1] + t**3 * pontos_controle[3][1]
                pontos_curva.append((x, y))
        else:
            # Para casos com mais pontos, usar método recursivo simplificado
            for i in range(num_segmentos + 1):
                t = i / num_segmentos
                # Algoritmo de De Casteljau simplificado
                pontos_temp = pontos_controle.copy()
                n = len(pontos_temp)
                for j in range(1, n):
                    for k in range(n - j):
                        pontos_temp[k] = (
                            (1 - t) * pontos_temp[k][0] + t * pontos_temp[k + 1][0],
                            (1 - t) * pontos_temp[k][1] + t * pontos_temp[k + 1][1]
                        )
                pontos_curva.append(pontos_temp[0])
        # Desenhar segmentos da curva usando Bresenham
        for i in range(len(pontos_curva) - 1):
            x0, y0 = pontos_curva[i]
            x1, y1 = pontos_curva[i+1]
            self.desenhar_linha_pixel(x0, y0, x1, y1)
        
        # Salva os pontos de controle para transformações
        self.vertices_2d_ativo = np.array(pontos_controle)
        self.desenhar_grid_no_canvas()

    # Funções de Solicitação de Entrada 
    
    def solicitar_reta(self):
        entrada = simpledialog.askstring("Desenhar Reta", "Insira (x1,y1,x2,y2):")
        if entrada:
            try:
                x1, y1, x2, y2 = map(int, entrada.split(','))
                self.desenhar_reta_bresenham(x1, y1, x2, y2)
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Use: x1,y1,x2,y2")

    def solicitar_circulo(self):
        entrada = simpledialog.askstring("Desenhar Círculo", "Insira o centro (x,y) e o raio:")
        if entrada:
            try:
                cx, cy, raio = map(int, entrada.split(','))
                self.desenhar_circulo(cx, cy, raio)
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Use: x,y,raio")

    def solicitar_bezier(self):
        entrada = simpledialog.askstring("Curva de Bézier", "Insira os pontos de controle (x1,y1,x2,y2,...):")
        if entrada:
            try:
                coords = list(map(float, entrada.split(',')))
                if len(coords) < 4 or len(coords) % 2 != 0:
                    raise ValueError("Número de coordenadas inválido")
                
                pontos = [(coords[i], coords[i+1]) for i in range(0, len(coords), 2)]
                if len(pontos) < 2:
                    raise ValueError("Precisa de pelo menos 2 pontos de controle")
                
                self.desenhar_bezier(pontos)
            except (ValueError, IndexError) as e:
                messagebox.showerror("Erro", f"Formato inválido: {e}. Use: x1,y1,x2,y2,...")

    def solicitar_polilinha(self):
        entrada = simpledialog.askstring("Polilinha", "Insira os pontos (x1,y1,x2,y2,...):")
        
        if entrada:
            try:
                coords = list(map(float, entrada.split(',')))
                vertices = np.array(coords).reshape(-1, 2)
                self.desenhar_poligono(vertices)
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Múltiplos de 2 coordenadas.")

    # Funções de Transformações
    def verificar_objeto(self):
        if self.vertices_ativos is None:
            messagebox.showwarning("Atenção", "Nenhum objeto para transformar. Desenhe uma forma primeiro.")
            return False
        return True

    def solicitar_transladar(self):
        if self.verificar_objeto():
            entrada = simpledialog.askstring("Translação", "Insira os valores de translação (dx,dy):")
            if entrada:
                try:
                    dx, dy = map(float, entrada.split(','))
                    matriz_transf = np.array([[1, 0], [0, 1]])
                    vetor_transf = np.array([dx, dy])
                    self.vertices_ativos = self.vertices_ativos @ matriz_transf + vetor_transf
                    self.redesenhar_forma()
                except (ValueError, IndexError):
                    messagebox.showerror("Erro", "Formato inválido. Use: dx,dy")

    def solicitar_rotacionar(self): 
        if self.verificar_objeto():
            entrada = simpledialog.askstring("Rotação", "Ângulo em graus e ponto de pivô (ângulo,px,py):")
            if entrada:
                try:
                    angulo, px, py = map(float, entrada.split(','))
            # 1. Transladar para a origem
                    vertices_transladados = self.vertices_ativos - np.array([px, py])
            # 2. Rotacionar em torno da origem
                    cos_a, sin_a = np.cos(angulo), np.sin(angulo)
                    matriz_rot = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
                    vertices_rotacionados = vertices_transladados @ matriz_rot
            # 3. Transladar de volta para a posição original
                    self.vertices_ativos = vertices_rotacionados + np.array([px, py])
                    self.redesenhar_forma()
                except (ValueError, IndexError):
                    messagebox.showerror("Erro", "Formato inválido. Use: ângulo,px,py")
  
    def solicitar_escalar(self):
        if self.verificar_objeto():
            entrada = simpledialog.askstring("Escala", "Fatores (sx,sy) e ponto fixo (sx,sy,px,py):")
            if entrada:
                try:
                    sx, sy, px, py = map(float, entrada.split(','))
                    matriz_esc = np.array([[sx, 0], [0, sy]])
                    self.vertices_ativos = self.vertices_ativos - np.array([px, py])
                    self.vertices_ativos = self.vertices_ativos @ matriz_esc
                    self.vertices_ativos = self.vertices_ativos + np.array([px, py])
                    self.vertices_ativos = self.vertices_ativos @ matriz_esc
                    self.redesenhar_forma()
                except (ValueError, IndexError):
                    messagebox.showerror("Erro", "Formato inválido. Use: sx,sy")

    # Funções de Processamentos

    def solicitar_preenchimento_recursivo(self):
       if len(self.vertices_ativos) == 1:
           entrada = simpledialog.askstring("Preenchimento", "Insira o ponto de partida (x,y):")
           if entrada:
                try:
                    x, y = map(float, entrada.split(','))
                    x_pixel, y_pixel = self.converter_para_canvas(x, y)
                    self.preencher(int(x_pixel), int(y_pixel), "red", "blue", "#f0f0f0")
                except (ValueError, IndexError):
                    messagebox.showerror("Erro", "Formato inválido. Use: x,y")

       else:
           self.limpar_formas()
           self.desenhar_poligono2(self.vertices_ativos)
           entrada = simpledialog.askstring("Preenchimento", "Insira o ponto de partida (x,y):")
           if entrada:
                try:
                    x, y = map(float, entrada.split(','))
                    x_pixel, y_pixel = self.converter_para_canvas(x, y)
                    self.preencher(int(x_pixel), int(y_pixel), "red", "blue", "#f0f0f0")
                except (ValueError, IndexError):
                    messagebox.showerror("Erro", "Formato inválido. Use: x,y")

    def preencher(self, x_pixel, y_pixel, cor_preenchimento, cor_borda, cor_fundo):
        """Preenchimento recursivo (Flood Fill) operando na matriz de pixels."""
        sys.setrecursionlimit(2000)

        fill_val = 2
        borda_val = 1
        fundo_val = 0
        
        stack = [(x_pixel, y_pixel)]

        while stack:
            px, py = stack.pop()
            
            if not (0 <= px < self.canvas_width and 0 <= py < self.canvas_height):
                continue
            if self.pixel_grid[px, py] != fundo_val:
                continue
            
            self.pixel_grid[px, py] = fill_val
            
            stack.append((px + 1, py))
            stack.append((px - 1, py))
            stack.append((px, py + 1))
            stack.append((px, py - 1))
        
        self.desenhar_grid_no_canvas()

    def solicitar_preenchimento_varredura(self):
        if len(self.vertices_ativos) == 1: # Se for Circulo
               borda_val = 1
               fill_val = 2
               # Percorre a grade de pixels da parte superior para a inferior e da esquerda para a direita
               for y in range(self.canvas_height):
                dentro = False
                for x in range(self.canvas_width):
                    if self.pixel_grid[x, y] == borda_val:
                        dentro = not dentro
                    # Se estiver dentro da forma (depois da primeira borda), preenche o pixel vermelho no canvas para exibir o processo
                    elif dentro:
                        self.pixel_grid[x, y] = fill_val
                        self.canvas.create_rectangle(x, y, x + 1, y + 1, fill="red", outline="red", tags=self.forma_ativa_tag)
            
                # Atualiza o canvas para exibir as mudanças em tempo real
                self.root.update_idletasks()
                self.root.update()
                
                self.desenhar_grid_no_canvas()
        # Ser for outras formas  
        else:               
           self.limpar_formas()
           self.desenhar_poligono2(self.vertices_ativos)

           # Repete o processo
           borda_val = 1
           fill_val = 2

           for y in range(self.canvas_height):
               dentro = False
               for x in range(self.canvas_width):
                    if self.pixel_grid[x, y] == borda_val:
                        dentro = not dentro
                    elif dentro:
                       self.pixel_grid[x, y] = fill_val
                       self.canvas.create_rectangle(x, y, x + 1, y + 1, fill="red", outline="red", tags=self.forma_ativa_tag)

               self.root.update_idletasks()
               self.root.update()
  
           self.desenhar_grid_no_canvas()

    def solicitar_recorte(self):
        if self.vertices_ativos is None:
            messagebox.showwarning("Atenção", "Nenhum objeto para recortar.")
            return
        entrada = simpledialog.askstring("Recorte", "Insira a janela de recorte (xmin,ymin,xmax,ymax):")
        if entrada:
            try:
                xmin, ymin, xmax, ymax = map(float, entrada.split(','))
                x1, y1 = self.vertices_ativos[0]
                x2, y2 = self.vertices_ativos[1]
                
                resultado = self.cohen_sutherland_clip(x1, y1, x2, y2, xmin, ymin, xmax, ymax)
                
                if resultado:
                    self.limpar_canvas()
                    self.desenhar_reta_bresenham(resultado[0], resultado[1], resultado[2], resultado[3])
                else:
                    messagebox.showinfo("Recorte", "A linha está totalmente fora da janela.")
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Use: xmin,ymin,xmax,ymax")

    def cohen_sutherland_clip(self, x1, y1, x2, y2, xmin, ymin, xmax, ymax):
        """Implementa o algoritmo de Cohen-Sutherland para recorte de linhas."""
        # Códigos de região
        INSIDE = 0  # 0000
        LEFT = 1    # 0001
        RIGHT = 2   # 0010
        BOTTOM = 4  # 0100
        TOP = 8     # 1000
        
        # Função para calcular o código de região de um ponto
        def calcular_codigo(x, y):
            codigo = INSIDE
            if x < xmin:
                codigo |= LEFT
            elif x > xmax:
                codigo |= RIGHT
            if y < ymin:
                codigo |= BOTTOM
            elif y > ymax:
                codigo |= TOP
            return codigo
        
        # Calcular códigos iniciais
        codigo1 = calcular_codigo(x1, y1)
        codigo2 = calcular_codigo(x2, y2)
        
        aceito = False
        
        while True:
            # Se ambos os pontos estão dentro da janela
            if codigo1 == 0 and codigo2 == 0:
                aceito = True
                break
            # Se ambos os pontos estão fora da janela e na mesma região
            elif codigo1 & codigo2 != 0:
                break
            # Pelo menos um ponto está fora da janela
            else:
                # Pega o ponto que está fora da janela
                codigo_fora = codigo1 if codigo1 != 0 else codigo2
                
                # Encontrar interseção
                if codigo_fora & TOP:  # Ponto acima da janela
                    x = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1)
                    y = ymax
                elif codigo_fora & BOTTOM:  # Ponto abaixo da janela
                    x = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1)
                    y = ymin
                elif codigo_fora & RIGHT:  # Ponto à direita da janela
                    y = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1)
                    x = xmax
                elif codigo_fora & LEFT:  # Ponto à esquerda da janela
                    y = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1)
                    x = xmin
                
                # Substituir o ponto fora da janela
                if codigo_fora == codigo1:
                    x1, y1 = x, y
                    codigo1 = calcular_codigo(x1, y1)
                else:
                    x2, y2 = x, y
                    codigo2 = calcular_codigo(x2, y2)
        
        if aceito:
            return [x1, y1, x2, y2]
        else:
            return None
    
    def redesenhar_forma(self):
        # Se a forma for um círculo, redesenha de forma especial
        if len(self.vertices_ativos) == 1:
             self.desenhar_circulo(self.vertices_ativos[0, 0], self.vertices_ativos[0, 1], 5) # Raio fixo para o círculo transformado
        else:
            self.desenhar_poligono(self.vertices_ativos)

    def limpar_canvas(self):
        """Limpa o canvas e redefine o estado."""
        self.canvas.delete("all")
        self.vertices_ativos = None
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)
        self.desenhar_plano_cartesiano()

if __name__ == "__main__":
    # Certifique-se de ter a biblioteca numpy instalada:
    # pip install numpy
    root = tk.Tk()
    app = Pratica_CG(root)
    root.mainloop()