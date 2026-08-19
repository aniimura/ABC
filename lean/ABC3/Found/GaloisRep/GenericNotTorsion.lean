import ABC3.Found.GaloisRep.TranslateAutAll
import ABC3.Found.GaloisRep.TorsionAll

/-!
# Galois (G5) 第 125 ブロック —— **★★★★★★生成点は捩れ点ではない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★第 65-72 の資産がそのまま効いた

第 118 の `exists_mulByNHom`(`[n]` の引き戻し)は
`n • (生成点) ≠ 0` を仮定として受けていた。★これを閉じる。

★★鍵は**自前の `exists_divisor_root`**(捩れの有限性のために積んだもの):

    n • P = 0  ⟹  ∃ k, 1 ≤ k, k ∣ n, ΨSq_k(x_P) = 0

★★★これを**関数体上の生成点に当てる**と `ΨSq_k(coordX) = 0` になるが、
`ΨSq_k` は `F` 係数の**非零多項式**であり、`coordX` は `F` 上超越的(第 116)——矛盾。

## ★★★★足場の記録

| 出所 | 使ったもの |
|---|---|
| 自前(第 65-72) | `exists_divisor_root`——捩れの有限性のために積んだ補題 |
| 自前(第 114) | `nonsingular_coord`——生成点が曲線の非特異点 |
| 自前(第 116) | `coordX_transcendental` |
| mathlib | `ΨSq_ne_zero`・`map_ΨSq`・`Polynomial.eval_map` |

★**捩れの有限性のために積んだ補題が、そのまま生成点の非捩れ性に効いた。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `genericPoint_not_torsion` | ★★★★★★**生成点は捩れ点でない** |
| `genericPoint_not_torsion_charZero` | ★★標数 0 版 |
| `exists_mulByNHom_charZero` | ★★★★★★**`[n]` の引き戻しが仮定なしで存在** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★★★★**生成点は捩れ点ではない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`n • G = 0` なら `ΨSq_k(coordX) = 0`(自前の `exists_divisor_root`)。
★★`ΨSq_k` は `F` 係数の非零多項式なので、`coordX` の超越性(第 116)に矛盾する。 -/
theorem genericPoint_not_torsion [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ∣ n → (k : F) ≠ 0) :
    n • genericPoint W ≠ 0 := by
  intro hzero
  obtain ⟨k, hk1, hkn, hkroot⟩ :=
    exists_divisor_root (W.map (algebraMap F W.FunctionField)) (nonsingular_coord W) n hn hzero
  rw [WeierstrassCurve.map_ΨSq, Polynomial.eval_map] at hkroot
  exact coordX_transcendental W (W.ΨSq_ne_zero (by exact_mod_cast hchar k hk1 hkn)) hkroot

/-- ★★標数 0 版。 -/
theorem genericPoint_not_torsion_charZero [DecidableEq F] [CharZero F]
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n) :
    n • genericPoint W ≠ 0 :=
  genericPoint_not_torsion W n hn (fun k hk _ => Nat.cast_ne_zero.2 (by omega))

/-- ★★★★★★**`[n]` の引き戻しが仮定なしで存在する**(標数 0)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 118 の `exists_mulByNHom` の仮定が本ブロックで外れた。 -/
theorem exists_mulByNHom_charZero [DecidableEq F] [CharZero F]
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n) :
    ∃ (x y : W.FunctionField) (h : (W.map (algebraMap F W.FunctionField)).Nonsingular x y),
      n • genericPoint W = Point.some x y h ∧
      ∃ μ : W.CoordinateRing →+* W.FunctionField,
        μ (genX W) = x ∧ μ (genY W) = y :=
  exists_mulByNHom W n (genericPoint_not_torsion_charZero W n hn)

/-! ## ★出典の紐付け(`.src`) -/

def genericPoint_not_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——生成点が捩れ点でないこと)",
    sectionId := "genell-thm-3-8" }

def exists_mulByNHom_charZero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——[n] の引き戻しが仮定なしで存在すること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
