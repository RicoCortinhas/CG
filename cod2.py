import tkinter as tk
import numpy as np
import sys
from tkinter import simpledialog, messagebox

# problemas:
# escalar só para +;
# rotacionar sem o ponto de pivo;
# não tem o preenchimento e recorte;
# curva de bezier só para numeros fixos
#

class Pratica_CG:
    def __init__(self, root):
        self.root = root
        self.root.title("Trabalho de Computação Gráfica")
        
        # --- Configurações do Canvas e Escala ---
        self.canvas_width = 999
        self.canvas_height = 500
        self.canvas = tk.Canvas(root, width=self.canvas_width, height=self.canvas_height, bg="#f0f0f0")
        self.canvas.pack(pady=10)
        
        # Fator de escala para converter coordenadas do plano para pixels
        self.escala = 5
        self.center_x = self.canvas_width / 2
        self.center_y = self.canvas_height / 2

        # matriz de pixels que armazena o estado da imagem
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)
        
        self.desenhar_plano_cartesiano()
        
        # Variável para armazenar os vértices da forma 2D ativa
        self.vertices_2d_ativo = None
        self.forma_ativa_tag = "forma"

        # --- Frames para a organização dos botões ---
        self.frame_botoes_forma = tk.Frame(root)
        self.frame_botoes_forma.pack(pady=5)
        
        tk.Label(self.frame_botoes_forma, text="Desenhar Forma:").pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Reta (Bresenham)", command=self.solicitar_reta).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Círculo", command=self.solicitar_circulo).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Curva de Bézier", command=self.solicitar_bezier).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Polilinha", command=self.solicitar_polilinha).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_forma, text="Quadrado", command=self.solicitar_quadrado).pack(side="left", padx=5)

        self.frame_botoes_transf = tk.Frame(root)
        self.frame_botoes_transf.pack(pady=5)
        
        tk.Label(self.frame_botoes_transf, text="Transformações:").pack(side="left", padx=5)
        tk.Button(self.frame_botoes_transf, text="Translação", command=self.solicitar_transladar).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_transf, text="Rotação", command=self.solicitar_rotacionar).pack(side="left", padx=5)
        tk.Button(self.frame_botoes_transf, text="Escala", command=self.solicitar_escalar).pack(side="left", padx=5)

        self.frame_processamento = tk.Frame(root)
        self.frame_processamento.pack(pady=5)
        
        #tk.Label(self.frame_processamento, text="Processamento:").pack(side="left", padx=5)
        tk.Button(self.frame_processamento, text="Preenchimento Recursivo", command=self.solicitar_preenchimento_recursivo).pack(side="left", padx=5)
        tk.Button(self.frame_processamento, text="Preenchimento Varredura", command=self.solicitar_preenchimento_varredura).pack(side="left", padx=5)
        tk.Button(self.frame_processamento, text="Recorte", command=self.solicitar_recorte).pack(side="left", padx=5)
        
        self.btn_limpar = tk.Button(root, text="Limpar", command=self.limpar_canvas)
        self.btn_limpar.pack(pady=10)

    # --- Funções do Plano Cartesiano ---
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
    

    # --- Funções para Desenho de Formas ---
    def limpar_formas(self):
        """Limpa apenas as formas desenhadas, mantendo a grade."""
        self.canvas.delete(self.forma_ativa_tag)
        self.vertices_2d_ativo = None
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)

    def desenhar_reta_bresenham(self, x0, y0, x1, y1):
        """Implementa o algoritmo de Bresenham para desenhar uma linha."""
        #self.limpar_formas()
                
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
        
        self.vertices_2d_ativo = np.array([[x0, y0], [x1, y1]])
        self.desenhar_grid_no_canvas()
        
    def desenhar_poligono(self, vertices):
        """Desenha um polígono a partir de uma lista de vértices [x, y]."""
        #inicio do codigo antigo
        self.limpar_formas()
        # Converte as coordenadas do plano cartesiano para pixels do canvas
        pontos = []
        for x, y in vertices:
            pontos.append(self.center_x + x * self.escala)
            pontos.append(self.center_y - y * self.escala)
        self.canvas.create_polygon(pontos, outline="blue", fill="", width=2, tags=self.forma_ativa_tag)
        self.vertices_2d_ativo = vertices
        self.vertices_2d_ativo = np.array([[x0, y0], [x1, y1]]) #novo
        #fim do codigo antigo
        
    
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
        #self.canvas.create_oval(x1, y1, x2, y2, outline="blue", fill="", width=2, tags=self.forma_ativa_tag)
        # Salva o centro do círculo como vértice para transformações
        self.vertices_2d_ativo = np.array([[cx, cy]])

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


    def desenhar_bezier(self, pontos_controle):
        """Desenha uma curva de Bézier na matriz de pixels."""
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)
        num_segmentos = 2
        
        pontos_curva = []
        for i in range(num_segmentos + 1):
            t = i / num_segmentos
            p_x = (1 - t)**2 * pontos_controle[0][0] + 2 * (1 - t) * t * pontos_controle[1][0] + t**2 * pontos_controle[2][0]
            p_y = (1 - t)**2 * pontos_controle[0][1] + 2 * (1 - t) * t * pontos_controle[1][1] + t**2 * pontos_controle[2][1]
            pontos_curva.append((p_x, p_y))

        for i in range(len(pontos_curva) - 1):
            x0, y0 = pontos_curva[i]
            x1, y1 = pontos_curva[i+1]
            self.desenhar_reta_bresenham(x0, y0, x1, y1)
        
        self.desenhar_grid_no_canvas()
    

    # --- Funções de Solicitação de Entrada ---
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
        entrada = simpledialog.askstring("Curva de Bézier", "Insira 3 pontos (x1,y1,x2,y2,x3,y3):")
        if entrada:
            try:
                coords = list(map(float, entrada.split(',')))
                pontos = [(coords[0], coords[1]), (coords[2], coords[3]), (coords[4], coords[5])]
                self.desenhar_bezier(pontos)
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Use: x1,y1,x2,y2,x3,y3")

    def solicitar_polilinha(self):
        entrada = simpledialog.askstring("Polilinha", "Insira os pontos (x1,y1,x2,y2,...):")
        if entrada:
            try:
                coords = list(map(float, entrada.split(',')))
                vertices = np.array(coords).reshape(-1, 2)
                self.desenhar_poligono(vertices)
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Múltiplos de 2 coordenadas.")

    def solicitar_quadrado(self):
        entrada = simpledialog.askstring("Desenhar Quadrado", "Insira (x,y) do canto inferior esquerdo e o lado:")
        if entrada:
            try:
                x, y, lado = map(int, entrada.split(','))
                vertices = np.array([[x, y], [x + lado, y], [x + lado, y + lado], [x, y + lado]])
                self.desenhar_poligono(vertices)
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Use: x,y,lado")


    # Funções de Transformações
    def verificar_objeto(self):
        if self.vertices_2d_ativo is None:
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
                    self.vertices_2d_ativo = self.vertices_2d_ativo @ matriz_transf + vetor_transf
                    self.redesenhar_forma()
                except (ValueError, IndexError):
                    messagebox.showerror("Erro", "Formato inválido. Use: dx,dy")

    #def solicitar_rotacionar(self):
    #    if self.verificar_objeto():
    #        entrada = simpledialog.askstring("Rotação", "Insira o ângulo em graus:")
    #        if entrada:
    #            try:
    #                angulo = np.radians(float(entrada))
    #                cos_a, sin_a = np.cos(angulo), np.sin(angulo)
    #                matriz_rot = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    #                self.vertices_2d_ativo = self.vertices_2d_ativo @ matriz_rot
    #                self.redesenhar_forma()
    #            except ValueError:
    #                messagebox.showerror("Erro", "Ângulo inválido.")
    def solicitar_rotacionar(self):
        if self.verificar_objeto():
        # Pede o ângulo de rotação
            entrada_angulo = simpledialog.askstring("Rotação", "Insira o ângulo em graus:")
            if not entrada_angulo:
                return

        # Pede o ponto de rotação (x, y)
            entrada_ponto_x = simpledialog.askstring("Ponto de Rotação", "Insira a coordenada X do ponto de rotação:")
            entrada_ponto_y = simpledialog.askstring("Ponto de Rotação", "Insira a coordenada Y do ponto de rotação:")
        
            if not entrada_ponto_x or not entrada_ponto_y:
                return

            try:
                angulo = np.radians(float(entrada_angulo))
                ponto_rotacao_x = float(entrada_ponto_x)
                ponto_rotacao_y = float(entrada_ponto_y)

            # 1. Transladar para a origem
                vertices_transladados = self.vertices_2d_ativo - np.array([ponto_rotacao_x, ponto_rotacao_y])

            # 2. Rotacionar em torno da origem
                cos_a, sin_a = np.cos(angulo), np.sin(angulo)
                matriz_rot = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
                vertices_rotacionados = vertices_transladados @ matriz_rot

            # 3. Transladar de volta para a posição original
                self.vertices_2d_ativo = vertices_rotacionados + np.array([ponto_rotacao_x, ponto_rotacao_y])
            
                self.redesenhar_forma()
            
            except ValueError:
                messagebox.showerror("Erro", "Ângulo ou ponto de rotação inválido.")

    def solicitar_escalar(self):
        if self.verificar_objeto():
            entrada = simpledialog.askstring("Escala", "Insira os fatores (sx,sy):")
            if entrada:
                try:
                    sx, sy = map(float, entrada.split(','))
                    matriz_esc = np.array([[sx, 0], [0, sy]])
                    self.vertices_2d_ativo = self.vertices_2d_ativo @ matriz_esc
                    self.redesenhar_forma()
                except (ValueError, IndexError):
                    messagebox.showerror("Erro", "Formato inválido. Use: sx,sy")

    # Funções de Processamentos

    def solicitar_preenchimento_recursivo(self):
        entrada = simpledialog.askstring("Preenchimento", "Insira o ponto de partida (x,y):")
        if entrada:
            try:
                x, y = map(float, entrada.split(','))
                x_pixel, y_pixel = self.converter_para_canvas(x, y)
                self.preencher_recursivo(int(x_pixel), int(y_pixel), "red", "blue", "#f0f0f0")
            except (ValueError, IndexError):
                messagebox.showerror("Erro", "Formato inválido. Use: x,y")

    def solicitar_preenchimento_varredura(self):
        entrada = simpledialog.askstring("Preenchimento", "Insira o ponto de partida (x,y):")

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



    # Funções de Execuções

    def preencher_recursivo(self, x_pixel, y_pixel, cor_preenchimento, cor_borda, cor_fundo):
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
    
    def redesenhar_forma(self):
        """Redesenha a forma ativa após uma transformação."""
        # Se a forma for um círculo, redesenha de forma especial
        if len(self.vertices_2d_ativo) == 1:
             self.desenhar_circulo(self.vertices_2d_ativo[0, 0], self.vertices_2d_ativo[0, 1], 5) # Raio fixo para o círculo transformado
        else:
            self.desenhar_poligono(self.vertices_2d_ativo)
    
    def limpar_canvas(self):
        """Limpa o canvas e redefine o estado."""
        self.canvas.delete("all")
        self.vertices_2d_ativo = None
        self.pixel_grid = np.zeros((self.canvas_width, self.canvas_height), dtype=int)
        self.desenhar_plano_cartesiano()

if __name__ == "__main__":
    # Certifique-se de ter a biblioteca numpy instalada:
    # pip install numpy
    root = tk.Tk()
    app = Pratica_CG(root)
    root.mainloop()