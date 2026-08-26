import ABC3.Found.GaloisRep.GenericNotTorsion
import ABC3.Found.GaloisRep.MulPoint

/-!
# Galois (G5) 第 323 ブロック —— **★★★★★★★生成点での乗法公式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★到達点

> **`x([n]P)·ΨSq_n(x) = Φ_n(x)`**——関数体 `F(E)` の生成点で(`mulByN_coordX`)

★★★これが第 196-198 で「非退化性に残る唯一の入力」と測定されたものである。

## ★★★★★在庫だけで組み上がった

| 段 | 在庫 |
|---|---|
| 点の水準の乗法公式 `x'·ΨSq_n(x) = Φ_n(x)` | ★第 52 `MulPoint.lean` の `mulOK_of_ne` |
| 仮定 `ΨSq_k(x) ≠ 0` (`1 ≤ k ≤ n`) | ★★第 116 の `coordX_transcendental` + mathlib の `ΨSq_ne_zero` |
| 生成点が非特異点であること | ★第 114 の `nonsingular_coord` |
| `[n]` の環準同型 `μ` | ★第 118 の `pointHom` |

★★★★**新しい数学は要らなかった**——`mulOK_of_ne` を**生成点に当てる**だけである。
★`ΨSq_k(coordX) ≠ 0` は「生成点は捩れ点でない」(第 125)と同じ 1 行で出る:
`ΨSq_k` は `F` 係数の非零多項式、`coordX` は `F` 上超越的。

## ★★★★★★これが (G5) のどこに効くか

    (G5) 非退化性
      → hfix : F(E)^{E[n]} = [n]^* F(E)          (第 197)
      → deg[n] = n²                              (第 196 の Artin が下から)
      → [F(x) : F(x∘[n])] <= n²                  (第 198:Φ_n − c·ΨSq_n はモニックで次数 n²)
      → **Φ_n(x) = x_n·ΨSq_n(x)**                ★本ブロック

★★残るのは `x_n` の超越性(`μ` の単射性に要る)と、体の拡大次数の帳簿だけである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psiSq_coordX_ne_zero` | ★★`ΨSq_k(coordX) ≠ 0`(`1 ≤ k`、標数 0) |
| `mulByN_coordX` | ★★★★★★**`x_n·ΨSq_n(x) = Φ_n(x)`** |
| `mulByN_coordX_hom` | ★★★★★★★環準同型 `μ` つきの形 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★分点多項式は生成点で消えない -/

/-- ★★**`ΨSq_k(coordX) ≠ 0`**——`ΨSq_k` は `F` 係数の非零多項式、`coordX` は超越的。

★第 125 の「生成点は捩れ点でない」と同じ 1 行である。 -/
theorem psiSq_coordX_ne_zero (W : WeierstrassCurve.Affine F) [CharZero F] (k : ℕ) (hk : 1 ≤ k) :
    Polynomial.eval (coordX W) ((W.map (algebraMap F W.FunctionField)).ΨSq (k : ℤ)) ≠ 0 := by
  rw [WeierstrassCurve.map_ΨSq, Polynomial.eval_map]
  refine coordX_transcendental W (W.ΨSq_ne_zero ?_)
  rw [Int.cast_natCast]
  exact Nat.cast_ne_zero.2 (by omega)

/-! ## ★★★★★★生成点での乗法公式 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**生成点での乗法公式** `x([n]P)·ΨSq_n(x) = Φ_n(x)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 52 の `mulOK_of_ne` を**生成点に当てる**だけで出る。 -/
theorem mulByN_coordX (W : WeierstrassCurve.Affine F) [W.IsElliptic] [DecidableEq F]
    [CharZero F] (n : ℕ) (hn : 1 ≤ n) :
    ∃ (x' y' : W.FunctionField) (h' : (W.map (algebraMap F W.FunctionField)).Nonsingular x' y'),
      n • genericPoint W = Point.some x' y' h'
        ∧ x' * Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ))
            = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.Φ (n : ℤ)) := by
  classical
  obtain ⟨x', y', h', heq, hx, _⟩ :=
    mulOK_of_ne (W.map (algebraMap F W.FunctionField)) (nonsingular_coord W) n hn
      (fun k hk1 _ => psiSq_coordX_ne_zero W k hk1)
  refine ⟨x', y', h', heq, ?_⟩
  rw [WeierstrassCurve.map_ΨSq, Polynomial.eval_map, WeierstrassCurve.map_Φ,
    Polynomial.eval_map] at hx
  exact hx

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**環準同型 `μ` つきの形**——`μ (genX W)` がちょうど `x([n]P)` である。

★`μ` は第 118 の `pointHom` である。 -/
theorem mulByN_coordX_hom (W : WeierstrassCurve.Affine F) [W.IsElliptic] [DecidableEq F]
    [CharZero F] (n : ℕ) (hn : 1 ≤ n) :
    ∃ (x' y' : W.FunctionField) (h' : (W.map (algebraMap F W.FunctionField)).Nonsingular x' y'),
      n • genericPoint W = Point.some x' y' h'
        ∧ x' * Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ))
            = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.Φ (n : ℤ))
        ∧ ∃ μ : W.CoordinateRing →+* W.FunctionField,
            μ (genX W) = x' ∧ μ (genY W) = y' := by
  classical
  obtain ⟨x', y', h', heq, hx, _⟩ :=
    mulOK_of_ne (W.map (algebraMap F W.FunctionField)) (nonsingular_coord W) n hn
      (fun k hk1 _ => psiSq_coordX_ne_zero W k hk1)
  rw [WeierstrassCurve.map_ΨSq, Polynomial.eval_map, WeierstrassCurve.map_Φ,
    Polynomial.eval_map] at hx
  exact ⟨x', y', h', heq, hx, pointHom W h'.1, pointHom_genX W h'.1, pointHom_genY W h'.1⟩

/-! ## ★出典の紐付け(`.src`) -/

def mulByN_coordX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——乗法 [n] の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

def mulByN_coordX_hom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——乗法 [n] の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
