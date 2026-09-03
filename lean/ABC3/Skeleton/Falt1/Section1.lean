import ABC3.Meta.Claim
import ABC3.Interface.Falt1.Ramification

/-!
# [Falt1] Chapter I §1 —— Ramification theory for discrete valuation rings(2 項目)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)。物理 p.4-5(印字 p.257-258)。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。★§1 は 2/2。
-/

namespace ABC3.Skeleton.Falt1

open ABC3.Interface.Falt1

/-- **`[Falt1] Lemma 1.1`**。

内容 (Falt1 p.4、260dpi 目視。OCR 層は数式記号を激しく壊す(Ω→"K"、
⊗→"?"、⊂→"c" 等)ため地の文で写す——Tate・BK と同じ運用): For any
extension V ⊂ W, as above, the natural map Ω_V ⊗_V W → Ω_W is
injective, and its cokernel Ω_{W/V} has the same length as W/p^δW. -/
def lemma_1_1 (E : RamificationSetup) :=
  And.intro E.lem11_injective E.lem11_length_eq

def lemma_1_1.nonvacuous := lemma_1_1 RamificationSetup.example

def lemma_1_1.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 4, item := "Lemma 1.1", sectionId := "falt1-lemma-1-1" }

/-- **`[Falt1] Theorem 1.2`**。

内容 (Falt1 p.5、260dpi 目視。OCR 層は数式記号を激しく壊すため地の文で
写す——Interface/Falt1/Ramification.lean 冒頭の注記参照): If V ⊂ W is
any extension and W_n denotes the normalization of the tensor product
V_n ⊗_V W, then the differents δ(W_n/V_n) of W_n over V_n (or more
precisely of each factor of the semilocal ring W_n) converge to 0 for
n → ∞. -/
def theorem_1_2 (E : RamificationSetup) := E.thm12

def theorem_1_2.nonvacuous := theorem_1_2 RamificationSetup.example

def theorem_1_2.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 5, item := "Theorem 1.2", sectionId := "falt1-thm-1-2" }

end ABC3.Skeleton.Falt1
