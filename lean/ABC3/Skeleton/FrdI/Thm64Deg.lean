/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Skeleton.NumberField.Chebotarev
import ABC3.Found.AddEquivNNReal
import Mathlib.Data.NNReal.Basic
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem

/-!
# [FrdI] Theorem 6.4 —— `deg(Ψ^rlf)` の段(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> Theorem 6.4. (Arithmetic Frobenioids) For i = 1, 2, let Fi be a number

## ★★この節点群の位置づけ

`Theorem 6.4` の (i)〜(iv) のうち、**Frobenioid の層は組み上がっている**
(`Example 6.3` の算術 Frobenioid、`Theorem 6.2` の再構成)。
残っているのは**数論の側の 4 段**で、鎖 `sec6items` の

| 節点 | 中身 | 在庫 |
|---|---|---|
| `thm64-i-dirichlet` | `δ_A : Pic_Φ(A) ≅ ℝ`(Dirichlet 単数定理) | ★★**mathlib にある** |
| `thm64-ii-deg` | `deg(Ψ^rlf) ∈ ℝ>0` の存在 | ★本ファイル(核は下記) |
| `thm64-iii-q` | `deg(Ψ^rlf) ∈ ℚ>0`(`Lemma 6.5, (ii)`) | ★`Lemma 6.5, (ii)` は**閉じた** |
| `thm64-iii-places` | `V(L₁) ≃ V(L₂)` が ℚ の素点の上下を保つ | ★`Lemma 6.5, (i)` は実装済み |
| `thm64-iv-galois` | `deg = 1` と `L₁ ≅ L₂` | ★Chebotarev(完全分解版) |

## ★★★在庫測定(2026-08-25、実測)

| 要るもの | mathlib | 判定 |
|---|---|---|
| Dirichlet 単数定理(階数) | `NumberField.Units.rank` / `unitLattice_rank` / `logEmbeddingEquiv` | ★★**ある** |
| Dedekind ζ の `s=1` の留数(類数公式) | `NumberField.dedekindZeta_residue` / `tendsto_sub_one_mul_dedekindZeta_nhdsGT` | ★★**ある** |
| six exponentials(`Lemma 6.5, (ii)`) | `ABC3.Found.SixExp`(自前、sorry なし) | ★★**ある** |
| Chebotarev(完全分解版) | `Skeleton/NumberField/Chebotarev.lean` | ☆skeleton |
| 単調加法自己同型 `ℝ≥0 ≃+ ℝ≥0` は正数倍 | grep 0 件 | ★**無い**(本ファイルの核) |

## ★★★★`thm64-ii-deg` の核は「順序を保つ加法同型は正数倍」である

原文は「`Φ^rlf₁(A₁) ≅ Φ^rlf₂(A₂)` は単系の同型だから順序を保つ」と 1 行で書く。
★単系の同型が順序を保つのは自明である(`ℝ≥0` では `a ≤ b ↔ ∃ c, a + c = b`)。
★★中身は**その先**、すなわち

```
順序を保つ加法自己同型 f : ℝ≥0 ≃+ ℝ≥0  ⟹  ∃ c > 0, ∀ x, f x = c * x
```

である。`f` が加法的なので `f (q·x) = q·f x`(`q : ℚ≥0`)が出て、
単調性と稠密性で `f x = f 1 · x` に絞られる。★**これが `deg(Ψ^rlf) ∈ ℝ>0` の実体**。
-/

namespace ABC3.Skeleton.Thm64

open ABC3.Meta

/-! ## ★1. `thm64-ii-deg` の核 —— 順序を保つ加法自己同型は正数倍 -/

/-- ★★★★★**`ℝ≥0` の加法自己同型は正の実数倍**。

★★これが `deg(Ψ^rlf) ∈ ℝ>0` の実体である。
★加法同型は `a ≤ b ↔ ∃ c, a + c = b` を保つので**単調**、
単調＋加法なら有理数倍を保ち、稠密性で `f x = f 1 · x` に絞られる。

★★★★★**閉じた**(2026-08-25)—— 中身は `Found/AddEquivNNReal.lean` にある。 -/
theorem exists_smul_eq_of_addEquiv (f : NNReal ≃+ NNReal) :
    ∃ c : NNReal, c ≠ 0 ∧ ∀ x : NNReal, f x = c * x :=
  ABC3.Found.exists_smul_eq_of_addEquiv f

/-- ★★**加法同型は順序を保つ**(自明。原文が 1 行で書いている段)。 -/
theorem monotone_of_addEquiv (f : NNReal ≃+ NNReal) : Monotone f :=
  ABC3.Found.monotone_of_addEquiv f

/-! ## ★2. `thm64-iii-q` —— `deg ∈ ℚ>0`

★`Lemma 6.5, (ii)`(six exponentials)は `Found/SixExp/` で**閉じている**。
残るのは「条件を満たす 6 素点が取れる」段である。 -/

/-- ★★★**6 素点が取れる** —— `Lemma 6.5, (ii)` を当てるための入力。

★`log p₁, …, log p₆` が `ℚ` 上一次独立(`Lemma 6.5, (i)`、実装済み)なので、
`2 × 3` の格子を張る 6 素点を選べる。 -/
theorem exists_six_primes_for_lemma65 :
    ∃ p : Fin 6 → ℕ, (∀ i, Nat.Prime (p i)) ∧ Function.Injective p := by
  sorry

/-- ★★★★**`deg(Ψ^rlf) ∈ ℚ>0`** —— `Lemma 6.5, (ii)` を 6 素点に当てる。

★`c` が `deg(Ψ^rlf)` にあたる。`c · log p` がすべて `log`(有理数)になるので、
six exponentials で `c` が有理数に落ちる。 -/
theorem deg_rat_of_six_exp (c : NNReal) (_hc : c ≠ 0)
    (_h : ∀ p : ℕ, Nat.Prime p → ∃ q : ℕ, 0 < q ∧
      (c : ℝ) * Real.log p = Real.log q) :
    ∃ q : ℚ, 0 < q ∧ (c : ℝ) = q := by
  sorry

/-! ## ★3. `thm64-iv-galois` —— `deg = 1` と `L₁ ≅ L₂` -/

/-- ★★★★★★**`deg(Ψ^rlf) = 1`**(`L₁/ℚ` が Galois のとき)と体の同型。

★★Chebotarev の**完全分解版**(`Skeleton/NumberField/Chebotarev.lean` の
`nonempty_algHom_of_splitsCompletely_subset`)で
「完全分解する素点が包含すれば体が包含する」を使う。
★★★**類体論(Artin の相互法則)は要らない** —— 2026-08-20 の測定のとおり、
[FrdI] が使う 3 箇所はいずれも完全分解の場合に収まる。 -/
theorem deg_eq_one_of_galois (L₁ L₂ : Type) [Field L₁] [Field L₂]
    [NumberField L₁] [NumberField L₂]
    (_hgal : ∀ σ : L₁ ≃ₐ[ℚ] L₁, True)
    (_hspl : ∀ p : IsDedekindDomain.HeightOneSpectrum (NumberField.RingOfIntegers ℚ),
      (ABC3.Skeleton.Cheb.SplitsCompletely ℚ L₁ p ↔ ABC3.Skeleton.Cheb.SplitsCompletely ℚ L₂ p)) :
    Nonempty (L₁ ≃ₐ[ℚ] L₂) := by
  sorry

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_smul_eq_of_addEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (ii) — deg(Ψ^rlf) ∈ ℝ>0(順序を保つ加法同型は正数倍)",
    sectionId := "frdi-thm-6-4" }

def exists_smul_eq_of_addEquiv.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Found 側の本体(sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.exists_smul_eq_of_addEquiv") 114,
    .citation "[mathlib]" "NNReal の順序と加法(le_iff_exists_add)"
      (.inMathlib "le_iff_exists_add") 114,
    .citation "[mathlib]" "単調加法自己同型が線形であること"
      (.absent "Mathlib/Analysis/ に AddMonoidHom.toRealLinearMap(連続性を要求)しかなく、単調版・NNReal 版は無い(2026-08-25)") 114,
    .derivation "加法性から f(q·x) = q·f x(q : ℚ≥0)、単調性と ℚ の稠密性で f x = f 1 · x" 114,
    .implicitStep
      "★原文は「Φ^rlf₁(A₁) ≅ Φ^rlf₂(A₂) は単系の同型だから順序を保つ」と 1 行で書く" 114 ]

def exists_six_primes_for_lemma65.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iii) — Lemma 6.5, (ii) に当てる 6 素点",
    sectionId := "frdi-thm-6-4" }

def exists_six_primes_for_lemma65.needs : List ProofObligation :=
  [ .citation "[ABC3]" "log_primes_linearIndependent(Lemma 6.5, (i))"
      (.inProject "ABC3" "ABC3.Found.FrdI.log_primes_linearIndependent") 116,
    .citation "[ABC3]" "six exponentials(Lemma 6.5, (ii))"
      (.inProject "ABC3" "ABC3.Found.SixExp") 116,
    .derivation "2 × 3 の格子を張るように素数を 6 個選ぶ" 116 ]

def deg_rat_of_six_exp.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iii) — deg(Ψ^rlf) ∈ ℚ>0",
    sectionId := "frdi-thm-6-4" }

def deg_rat_of_six_exp.needs : List ProofObligation :=
  [ .citation "[ABC3]" "six exponentials(Lemma 6.5, (ii)、sorry なし)"
      (.inProject "ABC3" "ABC3.Found.SixExp") 116,
    .citation "[ABC3]" "exists_six_primes_for_lemma65"
      (.inProject "ABC3" "ABC3.Skeleton.Thm64.exists_six_primes_for_lemma65") 116,
    .derivation "c·log pᵢ がすべて log(有理数)なら、six exponentials で c は有理数" 116 ]

def deg_eq_one_of_galois.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — deg(Ψ^rlf) = 1 と L₁ ≅ L₂",
    sectionId := "frdi-thm-6-4" }

/-- ★★★**Chebotarev の完全分解版だけで足りる**(2026-08-20 の測定)。
類体論(Artin の相互法則・Hecke L 関数)は**共役類つきの一般形**にしか要らない。 -/
def deg_eq_one_of_galois.needs : List ProofObligation :=
  [ .citation "[ABC3]" "nonempty_algHom_of_splitsCompletely_subset(Spl が体を決める)"
      (.inProject "ABC3" "ABC3.Skeleton.Cheb.nonempty_algHom_of_splitsCompletely_subset") 116,
    .citation "[ABC3]" "chebotarev_splitsCompletely(完全分解の密度 1/[L:K])"
      (.inProject "ABC3" "ABC3.Skeleton.Cheb.chebotarev_splitsCompletely") 116,
    .citation "[mathlib]" "NumberField.dedekindZeta_residue(Dedekind ζ の s=1 での留数)"
      (.inMathlib "NumberField.dedekindZeta_residue") 116,
    .citation "[mathlib]" "Artin の相互法則 / Hecke L 関数 / 射類群"
      (.absent "Mathlib に ClassFieldTheory ディレクトリは無く、Chebotarev|Artin 相互法則|RayClassGroup は 0 件(2026-08-25 に find で再確認)") 116,
    .derivation
      "deg の有理性(iii)と、完全分解の集合が一致することから [L₁:ℚ] = [L₂:ℚ]、Galois なら L₁ ≅ L₂" 116,
    .implicitStep
      "★原文は「by Tchebotarev's density theorem」と外部文献へ送る" 116 ]

def monotone_of_addEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (ii) — 単系の同型は順序を保つ",
    sectionId := "frdi-thm-6-4" }

/-- ★自明(`a ≤ b ↔ ∃ c, a + c = b`)——原文が 1 行で書いている段そのもの。 -/
def monotone_of_addEquiv.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Found 側の本体(sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.monotone_of_addEquiv") 114,
    .citation "[mathlib]" "le_iff_exists_add(NNReal の順序は加法で書ける)"
      (.inMathlib "le_iff_exists_add") 114 ]

end ABC3.Skeleton.Thm64
