import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Algebra.Field.ZMod
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import ABC3.Meta.Claim

/-!
# `𝔰𝔩₂` への随伴作用の既約性 —— [GenEll] Lemma 3.1, (iv) の部品

★**これは原典の項目の実装ではない。** `Lemma 3.1, (iv)` を証明するために要る道具であり、
対応する番号付き項目が原典に無いので **`.src` は付けない**。
位置づけは `ResearchPaper/genell-goal.md` の段 1・柱 2。

## ★なぜこれが要るか(2026-08-16 に (iv) の証明を解析して特定した)

原文 (GenEll p.14) の (iv) の証明は **[Serre], Chapter IV, §3.4, Lemma 3** を引く。
その [Serre] は `0_Source` に無いので**自分で証明する**ことになる。中身は

> `J ⊆ SL₂(ℤ_l)` が**閉**部分群で `SL₂(𝔽_l)` へ全射なら `J = SL₂(ℤ_l)`(l ≥ 5)

であり、標準的な道は**合同フィルトレーション** `K_n = ker(SL₂(ℤ_l) → SL₂(ℤ/l^n))` に沿う帰納法:

| 段 | 内容 | 状態 |
|---|---|---|
| 1 | `K_n / K_{n+1} ≅ 𝔰𝔩₂(𝔽_l)`、共役作用が随伴表現に落ちる | **未** |
| 2 | `J ∩ K_n` の像が部分加群になる | **未** |
| 3 | ★**`𝔰𝔩₂(𝔽_l)` が既約**——部分加群は `0` か全体 | ★**本ファイル** |
| 4 | `0` でないことを `x^l` の類の計算で見る(★ここで `l ≥ 5` が `C(l,3) ≡ 0` の形で効く) | **未** |

## ★実装上の判断

mathlib の Lie 代数(`Mathlib/Algebra/Lie/Sl2.lean`)や表現論を経由せず、
**2×2 行列の座標計算**で書いた。理由:

- 標数 `l` の modular representation theory は mathlib に実質無い
  (探索範囲: `Mathlib` 全体を `modular representation` / `Representation.*char` で grep、**1 件**)
- 必要なのは `Ad(upper t)` と `Ad(lower t)` による不変性だけで、それは**行列の積 3 つ**で書ける
- ★`u(t)⁻¹ = u(−t)` なので**逆行列を作らずに済む**

## ★★仮定は `(2 : K) ≠ 0` だけである(証明して分かった)

原文の `l ≥ 5` は `Lemma 3.1` 全体に掛かる仮定だが、
**既約性そのものに要るのは `2 ≠ 0` だけ**であった——
`e` を取り出す段の係数が `2c`、`h` と `f` を分ける段の係数が `2` である。
★`l ≥ 5` が効くのは (iv) の別の段(上の 4)であって、ここではない。

## ★`Submodule` で述べる理由(最初 `AddSubgroup` で書いて失敗した)

`AddSubgroup` は**スカラー倍で閉じない**ので、`(−2c) • e ∈ V` から `e ∈ V` が出ない。
表現の部分対象は部分加群なので `Submodule K` で述べる。
★素体 `𝔽_l` 上では加法部分群は自動的に部分空間なので、(iv) への適用では何も失わない。
-/

namespace ABC3.Found.GenEll

open Matrix

section Adjoint

variable {K : Type*} [Field K]

/-- `𝔰𝔩₂` の基底 `e = !![0,1;0,0]`。 -/
def slE : Matrix (Fin 2) (Fin 2) K := !![0, 1; 0, 0]

/-- `𝔰𝔩₂` の基底 `f = !![0,0;1,0]`。 -/
def slF : Matrix (Fin 2) (Fin 2) K := !![0, 0; 1, 0]

/-- `𝔰𝔩₂` の基底 `h = !![1,0;0,-1]`。 -/
def slH : Matrix (Fin 2) (Fin 2) K := !![1, 0; 0, -1]

/-- 上三角の基本行列による随伴作用 `X ↦ u(t) · X · u(t)⁻¹`(`u(t)⁻¹ = u(−t)`)。 -/
def adU (t : K) (X : Matrix (Fin 2) (Fin 2) K) : Matrix (Fin 2) (Fin 2) K :=
  !![1, t; 0, 1] * X * !![1, -t; 0, 1]

/-- 下三角の基本行列による随伴作用。 -/
def adL (t : K) (X : Matrix (Fin 2) (Fin 2) K) : Matrix (Fin 2) (Fin 2) K :=
  !![1, 0; t, 1] * X * !![1, 0; -t, 1]

/-- ★跡 0 の行列に対する `adU` の明示形。 -/
theorem adU_traceZero (t a b c : K) :
    adU t !![a, b; c, -a] = !![a + t * c, b - 2 * a * t - c * t ^ 2; c, -(a + t * c)] := by
  ext i j
  fin_cases i <;> fin_cases j <;> (simp [adU]; try ring)

/-- ★跡 0 の行列に対する `adL` の明示形。 -/
theorem adL_traceZero (t a b c : K) :
    adL t !![a, b; c, -a] = !![a - t * b, b; c + 2 * a * t - b * t ^ 2, -(a - t * b)] := by
  ext i j
  fin_cases i <;> fin_cases j <;> (simp [adL]; try ring)

/-- 跡 0 の行列は `!![a, b; c, -a]` の形に書ける。 -/
theorem eq_traceZero_form (X : Matrix (Fin 2) (Fin 2) K) (h : X 0 0 + X 1 1 = 0) :
    X = !![X 0 0, X 0 1; X 1 0, -(X 0 0)] := by
  ext i j
  fin_cases i <;> fin_cases j <;> (simp; try linear_combination h)

end Adjoint

section Irreducible

variable {K : Type*} [Field K]

variable (V : Submodule K (Matrix (Fin 2) (Fin 2) K))
  (hU : ∀ (t : K) {X}, X ∈ V → adU t X ∈ V)
  (hL : ∀ (t : K) {X}, X ∈ V → adL t X ∈ V)
  (htr : ∀ X ∈ V, X 0 0 + X 1 1 = 0)

include hU htr in
/-- ★左下成分が `0` でない元が `V` にあれば `e ∈ V`。

`adU 1` との差を **2 回**取ると、左下成分が消え、右上だけが残る。 -/
theorem slE_mem_of_lowerLeft_ne (hK : (2 : K) ≠ 0)
    {X : Matrix (Fin 2) (Fin 2) K} (hX : X ∈ V) (hc : X 1 0 ≠ 0) :
    (slE : Matrix (Fin 2) (Fin 2) K) ∈ V := by
  set a := X 0 0 with ha
  set b := X 0 1 with hb
  set c := X 1 0 with hc'
  have hform : X = !![a, b; c, -a] := eq_traceZero_form X (htr X hX)
  -- 1 回目: `Y = adU 1 X - X = !![c, -2a-c; 0, -c]`
  have hY : adU 1 X - X ∈ V := V.sub_mem (hU 1 hX) hX
  have hYval : adU 1 X - X = !![c, -(2 * a) - c; 0, -c] := by
    rw [hform, adU_traceZero]
    ext i j
    fin_cases i <;> fin_cases j <;> (simp; try ring)
  -- 2 回目: `Z = adU 1 Y - Y = (-(2*c)) • e`
  rw [hYval] at hY
  have hZ : adU 1 !![c, -(2 * a) - c; 0, -c] - !![c, -(2 * a) - c; 0, -c] ∈ V :=
    V.sub_mem (hU 1 hY) hY
  have hZval : adU 1 !![c, -(2 * a) - c; 0, -c] - !![c, -(2 * a) - c; 0, -c]
      = (-(2 * c)) • (slE : Matrix (Fin 2) (Fin 2) K) := by
    have : (!![c, -(2 * a) - c; 0, -c] : Matrix (Fin 2) (Fin 2) K)
        = !![c, -(2 * a) - c; 0, -c] := rfl
    rw [show (!![c, -(2 * a) - c; 0, -c] : Matrix (Fin 2) (Fin 2) K)
        = !![c, -(2 * a) - c; (0 : K), -c] from rfl]
    rw [adU_traceZero 1 c (-(2 * a) - c) 0]
    ext i j
    fin_cases i <;> fin_cases j <;> (simp [slE]; try ring)
  rw [hZval] at hZ
  have hne : -(2 * c) ≠ 0 := neg_ne_zero.2 (mul_ne_zero hK hc)
  have := V.smul_mem (-(2 * c))⁻¹ hZ
  rwa [smul_smul, inv_mul_cancel₀ hne, one_smul] at this

include hU hL htr in
/-- ★★**随伴表現の既約性**(基本行列による不変性だけを仮定した形)。

`V` が部分加群で、`adU t` と `adL t`(すべての `t`)で不変、跡 0 で、`0` でないなら、
`V` は `e`・`f`・`h` をすべて含む。

★仮定は `(2 : K) ≠ 0` だけ——`l ≥ 5` は要らない。 -/
theorem sl2_adjoint_irreducible (hK : (2 : K) ≠ 0)
    {X₀ : Matrix (Fin 2) (Fin 2) K} (hX₀ : X₀ ∈ V) (hX₀ne : X₀ ≠ 0) :
    (slE : Matrix (Fin 2) (Fin 2) K) ∈ V ∧ (slF : Matrix (Fin 2) (Fin 2) K) ∈ V
      ∧ (slH : Matrix (Fin 2) (Fin 2) K) ∈ V := by
  -- 段 1: 左下成分が `0` でない元を `V` の中に作る
  have hlow : ∃ Y ∈ V, Y 1 0 ≠ 0 := by
    by_cases hc : X₀ 1 0 = 0
    case neg => exact ⟨X₀, hX₀, hc⟩
    set a := X₀ 0 0 with ha
    set b := X₀ 0 1 with hb
    have hform : X₀ = !![a, b; 0, -a] := by
      rw [eq_traceZero_form X₀ (htr X₀ hX₀), hc]
    by_cases hb0 : b = 0
    · -- `X₀ = a·h`(`a ≠ 0`)。`adL 1` との差の左下は `2a ≠ 0`
      have ha0 : a ≠ 0 := by
        intro h
        exact hX₀ne (by rw [hform, h, hb0]; ext i j; fin_cases i <;> fin_cases j <;> simp)
      refine ⟨adL 1 X₀ - X₀, V.sub_mem (hL 1 hX₀) hX₀, ?_⟩
      rw [hform, adL_traceZero, hb0]
      simpa using mul_ne_zero hK ha0
    · -- `b ≠ 0`。左下は `t·(2a − b·t)` なので、`t ∈ {1, −1}` のどちらかで `0` でない
      have hex : ∃ t : K, t ≠ 0 ∧ 2 * a - b * t ≠ 0 := by
        by_cases h1 : 2 * a - b * 1 = 0
        · -- `2a = b` の場合。`t = −1` を使う
          refine ⟨-1, neg_ne_zero.2 one_ne_zero, ?_⟩
          intro h2
          -- `2a = b` かつ `2a = −b` なら `2b = 0`、よって `b = 0`
          apply hb0
          have e1 : 2 * a = b := by linear_combination h1
          have e2 : 2 * a = -b := by linear_combination h2
          have h4 : (2 : K) * b = 0 := by linear_combination e2 - e1
          rcases mul_eq_zero.1 h4 with h | h
          · exact absurd h hK
          · exact h
        · exact ⟨1, one_ne_zero, h1⟩
      obtain ⟨t, ht0, htne⟩ := hex
      refine ⟨adL t X₀ - X₀, V.sub_mem (hL t hX₀) hX₀, ?_⟩
      rw [hform, adL_traceZero]
      simp only [Matrix.sub_apply, Matrix.cons_val', Matrix.cons_val_zero, Matrix.empty_val',
        Matrix.cons_val_fin_one, Matrix.cons_val_one, Matrix.of_apply]
      have : (0 : K) + 2 * a * t - b * t ^ 2 - 0 = t * (2 * a - b * t) := by ring
      rw [this]
      exact mul_ne_zero ht0 htne
  obtain ⟨Y, hY, hYne⟩ := hlow
  have he : (slE : Matrix (Fin 2) (Fin 2) K) ∈ V :=
    slE_mem_of_lowerLeft_ne V hU htr hK hY hYne
  -- 段 2: `e` から `f` と `h` を作る
  -- `adL t e − e = !![-t, 0; -t², t]`
  have hA : adL 1 (slE : Matrix (Fin 2) (Fin 2) K) - slE ∈ V := V.sub_mem (hL 1 he) he
  have hB : adL (-1) (slE : Matrix (Fin 2) (Fin 2) K) - slE ∈ V := V.sub_mem (hL (-1) he) he
  have hAval : adL 1 (slE : Matrix (Fin 2) (Fin 2) K) - slE = !![-1, 0; -1, 1] := by
    show adL 1 (!![0, 1; 0, 0] : Matrix (Fin 2) (Fin 2) K) - !![0, 1; 0, 0] = _
    rw [show (!![0, 1; 0, 0] : Matrix (Fin 2) (Fin 2) K) = !![(0 : K), 1; 0, -(0 : K)] by
      ext i j; fin_cases i <;> fin_cases j <;> simp, adL_traceZero]
    ext i j
    fin_cases i <;> fin_cases j <;> (simp; try ring)
  have hBval : adL (-1) (slE : Matrix (Fin 2) (Fin 2) K) - slE = !![1, 0; -1, -1] := by
    show adL (-1) (!![0, 1; 0, 0] : Matrix (Fin 2) (Fin 2) K) - !![0, 1; 0, 0] = _
    rw [show (!![0, 1; 0, 0] : Matrix (Fin 2) (Fin 2) K) = !![(0 : K), 1; 0, -(0 : K)] by
      ext i j; fin_cases i <;> fin_cases j <;> simp, adL_traceZero]
    ext i j
    fin_cases i <;> fin_cases j <;> (simp; try ring)
  rw [hAval] at hA
  rw [hBval] at hB
  -- `A + B = (-2) • f`
  have hsum : (!![(-1 : K), 0; -1, 1] + !![(1 : K), 0; -1, -1])
      = (-2 : K) • (slF : Matrix (Fin 2) (Fin 2) K) := by
    ext i j
    fin_cases i <;> fin_cases j <;> (simp [slF]; try ring)
  have hf : (slF : Matrix (Fin 2) (Fin 2) K) ∈ V := by
    have h2 : (-2 : K) • (slF : Matrix (Fin 2) (Fin 2) K) ∈ V := by
      rw [← hsum]; exact V.add_mem hA hB
    have hne : (-2 : K) ≠ 0 := by
      simpa using hK
    have := V.smul_mem (-2 : K)⁻¹ h2
    rwa [smul_smul, inv_mul_cancel₀ hne, one_smul] at this
  -- `h = -A - f`
  have hh : (slH : Matrix (Fin 2) (Fin 2) K) ∈ V := by
    have : (slH : Matrix (Fin 2) (Fin 2) K)
        = -(!![(-1 : K), 0; -1, 1]) - (slF : Matrix (Fin 2) (Fin 2) K) := by
      ext i j
      fin_cases i <;> fin_cases j <;> (simp [slH, slF]; try ring)
    rw [this]
    exact V.sub_mem (V.neg_mem hA) hf
  exact ⟨he, hf, hh⟩

end Irreducible

end ABC3.Found.GenEll
