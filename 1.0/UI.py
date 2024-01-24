from PySide2.QtWidgets import *
import re
from sympy import Matrix, symbols
from sympy import solve_linear_system

app = QApplication([])

size = (670, 400)

main_window = QMainWindow()
main_window.resize(size[0], size[1])
main_window.move((1920 - size[0]) / 2, (1080 - size[1]) / 2)
main_window.setWindowTitle('Chem Equation Solver')

equation_text = QTextBrowser(main_window)
equation_text.append("Equation:")
equation_text.move(20, 20)
equation_text.resize(100, 30)

equation_input = QPlainTextEdit(main_window)
equation_input.setPlaceholderText("例（二硫键形成）：RSH -- (RS)2 + H[+] + e[-]")
equation_input.move(125, 20)
equation_input.resize(400, 30)

solve = QPushButton(main_window)
solve.setText("Solve")
solve.move(535, 20)
solve.resize(120, 30)

answer = QTextBrowser(main_window)
answer.move(20, 60)
answer.resize(440, 320)

tips = QTextBrowser(main_window)
tips.move(470, 60)
tips.resize(185, 320)
tips.setText("Wang125510 (C)\nDYC (C)\n因为没有使用多线程，\n在计算时无响应为正常情况，\n启动较慢问题是python的通病，\n本人目前仍未找到解决方法\n\n输入规范：\n下标数字直接写出，\n离子所带电荷写在方括号内跟随物质，\n物质间由\" + \"隔开，\n等号为\" -- \"，\n电子为\"e[-]\"")

al, acl, achl, anl, bl, bcl, bchl, bnl, elems, M = [], [], [], [], [], [], [], [], [], []


def error(e):
    answer.append("error " + e)
    # raise SystemExit


def char_check(equation):
    s = equation.strip().split(" -- ")
    if len(s) == 2:
        a, b = s[0], s[1]
    else:
        error("等号错误")

    re_str = "([qwertyuiopasdfghjklzxcvbnm\(\)QWERTYUIOPASDFGHJKLZXCVBNM\[\]1234567890+-]*\s*)*"
    res = re.findall(re_str, a)
    if len(res) != 2:
        error("反应物错误")
    res = re.findall(re_str, b)
    if len(res) != 2:
        error("生成物错误")
    global al, bl
    al, bl = a.split(" + "), b.split(" + ")
    return


def chem_check(chem):
    if chem == "e[-]":
        return
    re_str = "(\({0,1}[A-Z][a-z]{0,1}\d*\){0,1}){1,}(\[\d*[+-]\]){0,1}"
    res = re.findall(re_str, chem)
    if len(res) != 1:
        error("化学式错误")
    return


def spl_chem(chem):
    if chem == "e":
        return [("e", 1)]
    re_str = "[A-Z][a-z]{0,1}\d*"
    res = re.findall(re_str, chem)
    ls = []
    for i in res:
        if i[-1].isalpha():
            ls.append((i, 1))
        else:
            if i[1].isalpha():
                ls.append((i[0:2], int(i[2:])))
            else:
                ls.append((i[0], int(i[1:])))
    return ls


def spl(chem, n):
    if chem[-1] == "]":
        ls = chem[:-1].split("[")
        if len(ls[1]) == 1:
            try:
                charge = int(ls[1] + "1")
            except ValueError:
                error("电荷错误")
        else:
            try:
                charge = int(ls[1][-1] + ls[1][:-1])
            except ValueError:
                error("电荷错误")
        chem = ls[0]
    else:
        charge = 0
    if n == 0:
        global achl
        achl.append(charge)
    else:
        global bchl
        bchl.append(charge)

    bol = 0
    a, b, c = "", "", ""
    cl = []
    i = 0
    while i < len(chem):
        if chem[i] == "(":
            if bol == 0:
                bol = 1
                cl += spl_chem(a)
                a = ""
            else:
                error("括号匹配错误")
        elif chem[i] == ")":
            if bol == 1:
                bol = 0
                ls = spl_chem(b)
                b = ""
                j = i + 1
                c = ""
                m = 0
                while j < len(chem):
                    if chem[j].isdigit():
                        c += chem[j]
                        j += 1
                    else:
                        m = int(c)
                        i = j
                        break
                if m == 0:
                    m = int(c)
                for j in range(len(ls)):
                    ls[j] = (ls[j][0], ls[j][1] * m)
                cl += ls
            else:
                error("括号匹配错误")
        else:
            if bol == 0:
                a += chem[i]
            else:
                b += chem[i]
        i += 1
    if chem[-1] != ")":
        cl += spl_chem(a)

    cl_ = []
    bol = 0
    for i in range(len(cl)):
        for j in range(len(cl_)):
            if cl[i][0] == cl_[j][0]:
                cl_[j] = (cl_[j][0], cl_[j][1] + cl[i][1])
                bol = 1
                break
        if bol == 0:
            cl_.append(cl[i])

    if n == 0:
        global acl
        acl.append(cl_)
    else:
        global bcl
        bcl.append(cl_)
    return


def get_elem(cl):
    global elems
    for chem in cl:
        if chem != [("e", 1)]:
            for e in chem:
                if e[0] not in elems:
                    elems.append(e[0])
    return


def count_elem(acl, bcl):
    global M
    for e in elems:
        L = []
        for chem in acl:
            bol = 0
            for e_ in chem:
                if e_[0] == e:
                    L.append(e_[1])
                    bol = 1
                    break
            if bol == 0:
                L.append(0)
        for chem in bcl:
            bol = 0
            for e_ in chem:
                if e_[0] == e:
                    L.append(-e_[1])
                    bol = 1
                    break
            if bol == 0:
                L.append(0)
        L.append(0)
        M.append(L)
    return


def solve_lin_eqs(A):
    Mat = Matrix(A)
    count = Mat.shape[1] - 1
    sml = []
    for i in range(count):
        s = chr(i + 97)
        sml.append(symbols(s))
    cmd = "res = solve_linear_system(Mat, "
    for i in range(count):
        cmd += "sml[" + str(i) + "]" + ', '
    cmd = cmd[:-2] + ")"
    exec(cmd)
    rr = locals()['res']
    res_dict = {}
    for i, j in rr.items():
        res_dict[str(i)] = str(j)
    return res_dict


def cal(acl, achl, bcl, bchl):
    get_elem(acl)
    get_elem(bcl)
    count_elem(acl, bcl)
    L = achl
    for i in bchl:
        L.append(-i)
    L.append(0)
    M.append(L)
    res = solve_lin_eqs(M)
    answer.append(str(res))
    if all(value == "0" for value in res.values()):
        error("无解")
    deno = []
    for v in res.values():
        if '/' in v:
            deno.append(int(v.split('/')[-1]))
    if deno:
        num = deno.copy()
        while True:
            mi = min(num)
            ma = max(num)
            if mi != ma:
                ind = num.index(mi)
                num[ind] += deno[ind]
            else:
                break
        lcm = num[0]
    else:
        lcm = 1

    for alph in res['a']:
        if alph.isalpha():
            break
    if not alph.isalpha():
        answer.append(str(res))
        error("有无穷多解 或 存在一物质系数为0")
    ans = []
    for i, j in res.items():
        try:
            ans.append(int(eval(str(j).replace(alph, str(lcm)))))
        except NameError:
            answer.append(str(res))
            error("有无穷多解 或 存在一物质系数为0")
    ans.append(lcm)
    return ans


def main():
    global al, acl, achl, anl, bl, bcl, bchl, bnl, elems, M
    al, acl, achl, anl, bl, bcl, bchl, bnl, elems, M = [], [], [], [], [], [], [], [], [], []
    equation = equation_input.toPlainText().strip()
    char_check(equation)
    for i in al:
        chem_check(i)
        spl(i, 0)
    for i in bl:
        chem_check(i)
        spl(i, 1)
    ans = cal(acl, achl, bcl, bchl)
    nal, nbl = len(al), len(bl)
    equation = ""
    for i in range(nal):
        if ans[i] != 1:
            equation += str(ans[i])
        equation += al[i] + " + "
    equation = equation[:-3] + " == "
    for i in range(nbl):
        if ans[i + nal] != 1:
            equation += str(ans[i + nal])
        equation += bl[i] + " + "
    equation = equation[:-3]
    answer.append(equation)


solve.clicked.connect(main)

main_window.show()

app.exec_()
