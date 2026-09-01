/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38Bridge
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★`α` が像に入る段の行列表示（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★これは何か

`§9-992` で `Theorem 3.8` の最終段のうち**群論の 2 段**が閉じた。
残る 1 段は「`l ∤ v(q)` なら mod `l` 像が `α = (1 1 / 0 1)` を含む」である。

★★★**本ファイルはその段の「行列表示の側」を取る**——
Tate 一意化の言葉で `M_l(E) ⊆ Lˣ/q^ℤ` の元を `ζ^a·π^b`（`ζ^l = 1`、`π^l = q`）と書くと、

    `σ(ζ) = ζ`  かつ  `σ(π) = ζ·π`   ⟹   `σ(ζ^a π^b) = ζ^{a+b} π^b`

であり、座標 `(a,b)` は **`α = (1 1 / 0 1)` で動く**（`upperM_one_mulVec`）。

## ★残るのは Kummer の側（存在）

☆「`l ∤ v_K(q)` なら**そのような `σ` がある**」——`K(ζ_l, q^{1/l})/K(ζ_l)` の
Galois 群が `μ_l` と Kummer 双対であること。
★`[K(ζ_l):K] ∣ l−1` は `l` と素なので、`l ∤ v_K(q)` から `q ∉ K(ζ_l)^{×l}` が出る。

★★`lemma_3_2_i`（`Lemma32StableLine.lean`）は**逆向き**を取っている
——「安定な直線があって `l ∤ k` なら `l ∣ v_K(q)`」。
★★★本ファイルが要求するのは**その対偶ではなく、`σ` の存在**である。

## ★★★★★これで `Theorem 3.8` に残るのは

| 段 | 状態 |
|---|---|
| `α` が像に入る（行列表示） | ★★**本ファイル** |
| `α` が像に入る（Kummer の存在） | ☆残る |
| 安定な直線が無い ⟹ 非上三角 | ★`§9-992`（群論の側）／☆`l`-巡回 ⟷ 安定直線は残る |
| `Lemma 3.1, (iv)` | ★済み（`Sl2Padic.lean`） |
| `Lemma 3.7` | ☆`Prop 3.4` 待ち |
| `torsionExt` | ☆残る |

★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GenEll

open Matrix

/-! ## ★★★★★★★★★★★★★`σ` は座標を `α` で動かす -/

/-- ★★★★★★★★★★★★★**`σ(ζ) = ζ` かつ `σ(π) = ζπ` なら `σ(ζ^a π^b) = ζ^{a+b} π^b`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★Tate 一意化の言葉で `M_l(E) ⊆ Lˣ/q^ℤ` の元を `ζ^a·π^b` と書いたときの
`σ` の作用である（`ζ^l = 1`、`π^l = q`）。
★★座標 `(a,b)` は `(a+b, b)` へ動く——すなわち **`α = (1 1 / 0 1)`** である。 -/
theorem sigma_acts_as_alpha {K L : Type} [Field K] [Field L] [Algebra K L]
    (ζ π : Lˣ) (σ : L ≃ₐ[K] L)
    (hζ : Units.map (σ : L →* L) ζ = ζ)
    (hπ : Units.map (σ : L →* L) π = ζ * π) (a b : ℤ) :
    Units.map (σ : L →* L) (ζ ^ a * π ^ b) = ζ ^ (a + b) * π ^ b := by
  rw [map_mul, map_zpow, map_zpow, hζ, hπ, mul_zpow, zpow_add, mul_assoc]

/-- ★★★★★**座標の側** —— `α = (1 1 / 0 1)` は `(a,b) ↦ (a+b, b)` である。

★`sigma_acts_as_alpha` の指数の動きが、まさにこの行列の作用である。 -/
theorem upperM_one_mulVec {l : ℕ} (a b : ZMod l) :
    (upperM (1 : ZMod l)).mulVec ![a, b] = ![a + b, b] := by
  ext i
  fin_cases i <;> simp [upperM, Matrix.mulVec]

/-! ## ★★★★★★★★★★★★合同から行列の作用へ -/

/-- ★★★★★★★★★★★★★★
**座標の合同はそのまま `α` の作用**——★**無条件**（第 1175）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`a' ≡ a + b`、`b' ≡ b`（`mod l`）なら `ZMod l` の中で
`α · (a, b) = (a', b')` である。
★★★これが第 1174 の `tate_sigma_coord_alpha` の結論（`ℤ` の合同）を
**行列の言葉**に直す最後の 1 行である。 -/
theorem upperM_one_mulVec_of_dvd {l : ℕ} [NeZero l] {a b a' b' : ℤ}
    (ha : (l : ℤ) ∣ (a + b - a')) (hb : (l : ℤ) ∣ (b - b')) :
    (upperM (1 : ZMod l)).mulVec ![(a : ZMod l), (b : ZMod l)]
      = ![(a' : ZMod l), (b' : ZMod l)] := by
  rw [upperM_one_mulVec]
  have h1 : ((a : ZMod l) + (b : ZMod l)) = (a' : ZMod l) := by
    have h := (ZMod.intCast_zmod_eq_zero_iff_dvd (a + b - a') l).2 ha
    push_cast at h
    exact sub_eq_zero.mp h
  have h2 : (b : ZMod l) = (b' : ZMod l) := by
    have h := (ZMod.intCast_zmod_eq_zero_iff_dvd (b - b') l).2 hb
    push_cast at h
    exact sub_eq_zero.mp h
  rw [h1, h2]

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def upperM_one_mulVec_of_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標の合同はそのまま α の作用。★無条件)",
    sectionId := "genell-thm-3-8" }

def sigma_acts_as_alpha.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ(ζ)=ζ かつ σ(π)=ζπ なら座標は α で動く)",
    sectionId := "genell-thm-3-8" }

def upperM_one_mulVec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(α = (1 1 / 0 1) は (a,b) ↦ (a+b, b))",
    sectionId := "genell-thm-3-8" }

def sigma_acts_as_alpha.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_2_i(安定な直線の側——逆向き、第 410 ほか)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_2_i") 6,
    .folklore
      ("Kummer の存在: `l ∤ v_K(q)` なら σ(ζ) = ζ かつ σ(π) = ζ·π なる σ がある。" ++
       "★K(ζ_l, q^{1/l})/K(ζ_l) の Galois 群が μ_l と Kummer 双対であること。" ++
       "★★[K(ζ_l):K] ∣ l−1 は l と素なので l ∤ v_K(q) から q ∉ K(ζ_l)^{×l} が出る。**残る**") 6,
    .implicitStep
      ("★★★★★測定(2026-08-29): §9-992 で Theorem 3.8 の最終段のうち**群論の 2 段**が閉じ、" ++
       "本ファイルで**3 段目の行列表示の側**が取れた。" ++
       "★残るのは Kummer の**存在**の側だけである" ++
       "——lemma_3_2_i は逆向き(安定な直線があれば l ∣ v(q))を取っており、" ++
       "本件が要求するのは σ の存在である") 6,
    .implicitStep
      ("★★Theorem 3.8 の残高(2026-08-29): " ++
       "(1) α が像に入る Kummer の存在、(2) l-巡回 ⟷ 安定直線の対応、" ++
       "(3) Lemma 3.7(Prop 3.4 待ち)、(4) torsionExt(3・5 捩れ有理化 ⟹ 半安定)。" ++
       "★群論(Lemma 3.1, (iv)・§9-992・本ファイル)と Galois 表現の構成(galRep)は済んでいる") 7 ]

end ABC3.Found.GenEll
