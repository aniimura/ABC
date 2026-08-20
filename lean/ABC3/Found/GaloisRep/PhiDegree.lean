import ABC3.Found.GaloisRep.NondegStep

/-!
# Galois (G5) 第 198 ブロック —— **★★★★★★`Φ_n − c·ΨSq_n` はモニックで次数 `n²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★残る 1 行の形が定まった

第 197 で非退化性は **`F(E)^{E[n]} = [n]^* F(E)`** の 1 つに絞れた。
★第 196 の Artin から `[F(E) : F(E)^{E[n]}] = n²` なので、要るのは

    [F(x) : F(x∘[n])] ≤ n²

の**上からの評価だけ**である(下からは Artin が与える)。

### ★★★★★★上からの評価は「モニック多項式 1 本」で済む

`x` が `F(x_n)` 上で次数 `n²` のモニック多項式の根であれば良い。★その多項式は

    Φ_n(X) − x_n · ΨSq_n(X)

である。★★mathlib の `natDegree_Φ_le` と `coeff_Φ` から **`Φ_n` はモニックで次数ちょうど `n²`**
(仮定は何も要らない)。★★★`ΨSq_n` の次数は `n² − 1` 以下なので、差もモニックで次数 `n²`。

### ★★★★★★★これで「約分」の心配が消えた

当初は `finrank_eq_max_natDegree`(Lüroth の副産物)で `max(deg num, deg denom)` を使う
つもりだったが、それには **`Φ_n` と `ΨSq_n` が互いに素**であることが要った。
★モニック多項式による上からの評価に切り替えたので、**互いに素性は要らない**。
★★下からの評価は Artin が与えるので、上下が挟んで等号になる。

### ★残る入力(1 本だけ)

    Φ_n(x) = x_n · ΨSq_n(x)     (関数体 `F(E)` での等式)

★これが `Skeleton/GaloisRep/WeilFunctionField.lean` の `exists_mulByNPullback` である。
★★mathlib は分点多項式とその次数を持つが、**群法則とのリンクを持たない**
(`ΨSq`/`Φ` は `DivisionPolynomial/` の外に 0 件、2026-08-20 実測)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `natDegree_Φ_eq` | ★★★★**`deg Φ_n = n²` ちょうど**(仮定なし) |
| `monic_Φ` | ★★★★**`Φ_n` はモニック**(仮定なし) |
| `mulByNPoly` | `Φ_n − C c · ΨSq_n` |
| `monic_mulByNPoly` | ★★★★★★**モニックで次数 `n²`** |
-/

namespace ABC3.Found.GaloisRep

open Polynomial WeierstrassCurve

variable {F : Type} [Field F] (W : WeierstrassCurve F)

/-- ★★★★**`deg Φ_n = n²` ちょうど**——`natDegree_Φ_le` と `coeff_Φ` から。 -/
theorem natDegree_Φ_eq (n : ℤ) : (W.Φ n).natDegree = n.natAbs ^ 2 :=
  le_antisymm (W.natDegree_Φ_le n)
    (le_natDegree_of_ne_zero (by rw [W.coeff_Φ n]; exact one_ne_zero))

/-- ★★★★**`Φ_n` はモニック**。 -/
theorem monic_Φ (n : ℤ) : (W.Φ n).Monic := by
  rw [Polynomial.Monic, Polynomial.leadingCoeff, natDegree_Φ_eq W n, W.coeff_Φ n]

/-- `x` が `F(x_n)` 上で満たすべき多項式 `Φ_n(X) − c·ΨSq_n(X)`。 -/
noncomputable def mulByNPoly (n : ℤ) (c : F) : Polynomial F :=
  W.Φ n - Polynomial.C c * W.ΨSq n

/-- ★★★★★★**`Φ_n − c·ΨSq_n` はモニックで次数ちょうど `n²`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが `[F(x) : F(x∘[n])] ≤ n²` を与える多項式である。
★★`Φ_n` と `ΨSq_n` の互いに素性は**要らない**。 -/
theorem monic_mulByNPoly (n : ℤ) (hn : 1 ≤ n.natAbs) (c : F) :
    (mulByNPoly W n c).Monic ∧ (mulByNPoly W n c).natDegree = n.natAbs ^ 2 := by
  have hd := natDegree_Φ_eq W n
  have hmonic := monic_Φ W n
  have hsq : 1 ≤ n.natAbs ^ 2 := Nat.one_le_pow 2 _ (by omega)
  have hstep : (Polynomial.C c * W.ΨSq n).natDegree ≤ n.natAbs ^ 2 - 1 := by
    refine Polynomial.natDegree_mul_le.trans ?_
    rw [Polynomial.natDegree_C, zero_add]
    exact W.natDegree_ΨSq_le n
  have hlt : (Polynomial.C c * W.ΨSq n).natDegree < (W.Φ n).natDegree := by
    rw [hd]; omega
  exact ⟨hmonic.sub_of_left (Polynomial.degree_lt_degree hlt),
    (Polynomial.natDegree_sub_eq_left_of_natDegree_lt hlt).trans hd⟩

/-! ## ★出典の紐付け(`.src`) -/

def natDegree_Φ_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——Φ_n の次数)",
    sectionId := "genell-thm-3-8" }

def monic_mulByNPoly.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——deg[n] = n² のモニック多項式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
