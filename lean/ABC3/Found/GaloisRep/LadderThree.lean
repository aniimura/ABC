import ABC3.Found.GaloisRep.KummerPoint

/-!
# Galois (G5) 第 202 ブロック —— **★★★★★★梯子と EDS が噛み合うことを `n = 3` で確かめた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★梯子の 1 段は多項式恒等式になる

Montgomery 梯子の奇数段は `A = 2P`・`B = P`(差は `P`)で、第 201 の積の公式から

    x(3P)·x·(x(2P) − x)² = x(2P)²x² − b₄x(2P)x − b₆(x(2P)+x) − b₈

★`x(2P) = Φ₂/ΨSq₂`(第 199)を入れて分母を払うと、**`x` の多項式恒等式**になる。

### ★★★★★★鍵は `Φ₂ − X·ΨSq₂ = −Ψ₃`

`Φ_n := X·ΨSq_n − preΨ_{n+1}·preΨ_{n−1}·(…)` という定義から、`n = 2` では
`preΨ₃ = Ψ₃`・`preΨ₁ = 1` なので **`Φ₂ − X·ΨSq₂ = −Ψ₃`**。
★これで `(x(2P) − x)²` の分子が `Ψ₃²` になり、`ΨSq₃ = Ψ₃²` とちょうど打ち消す。

### ★★★★★残った恒等式

    Φ₃·X = Φ₂²X² − b₄Φ₂X·ΨSq₂ − b₆Φ₂·ΨSq₂ − b₆X·ΨSq₂² − b₈·ΨSq₂²

★`ring` で通った(2.6 秒)。★★**mathlib の `preΨ₄` の定義が群法則と噛み合っている**
ことの確認でもある——`Φ₃ = X·Ψ₃² − preΨ₄·Ψ₂Sq` の `preΨ₄` がここで効く。

## ★★この先

一般の `n` は次の 2 本の帰納法で回る:

| 段 | 使うもの |
|---|---|
| `n → 2n` | 2 倍公式(第 199) |
| `n → 2n+1` | Kummer の積の公式(第 201)と `x(P)` |

★どちらも「分母を払うと多項式恒等式」になり、その恒等式は mathlib の
`preΨ_even`・`preΨ_odd` から出る。★★`n = 3` で実際に通ったので、道は確実である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `Φ_two_sub_X_mul` | ★★★★★**`Φ₂ − X·ΨSq₂ = −Ψ₃`** |
| `ladder_three` | ★★★★★★**`n = 3` の梯子恒等式** |
-/

namespace ABC3.Found.GaloisRep

open Polynomial WeierstrassCurve

variable {F : Type} [Field F] (W : WeierstrassCurve F)

/-- ★★★★★**`Φ₂ − X·ΨSq₂ = −Ψ₃`**——`Φ` の定義そのもの(`preΨ₃ = Ψ₃`、`preΨ₁ = 1`)。 -/
theorem Φ_two_sub_X_mul : W.Φ 2 - Polynomial.X * W.ΨSq 2 = -W.Ψ₃ := by
  rw [WeierstrassCurve.Φ_two, WeierstrassCurve.ΨSq_two]
  simp only [WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.Ψ₃, WeierstrassCurve.b₂,
    WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈,
    map_ofNat, C_add, C_sub, C_mul, C_pow]
  ring

/-- ★★★★★★**`n = 3` の梯子恒等式**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 201 の Kummer の積の公式に `x(2P) = Φ₂/ΨSq₂` を入れて分母を払った形。
★★`Φ₂ − X·ΨSq₂ = −Ψ₃` と `ΨSq₃ = Ψ₃²` が打ち消し合う。 -/
theorem ladder_three :
    W.Φ 3 * Polynomial.X
      = (W.Φ 2) ^ 2 * Polynomial.X ^ 2 - Polynomial.C W.b₄ * (W.Φ 2) * Polynomial.X * (W.ΨSq 2)
        - Polynomial.C W.b₆ * (W.Φ 2) * (W.ΨSq 2)
        - Polynomial.C W.b₆ * Polynomial.X * (W.ΨSq 2) ^ 2
        - Polynomial.C W.b₈ * (W.ΨSq 2) ^ 2 := by
  rw [WeierstrassCurve.Φ_two, WeierstrassCurve.ΨSq_two, WeierstrassCurve.Φ_three]
  simp only [WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.Ψ₃, WeierstrassCurve.preΨ₄,
    WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈,
    map_ofNat, C_add, C_sub, C_mul, C_pow]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def Φ_two_sub_X_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——Φ₂ − X·ΨSq₂ = −Ψ₃)",
    sectionId := "genell-thm-3-8" }

def ladder_three.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——n = 3 の梯子恒等式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
