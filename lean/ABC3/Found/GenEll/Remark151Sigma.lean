/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.RatTowerLocalization
import ABC3.Found.GenEll.LogCondSigma

/-!
# [GenEll] Remark 1.5.1 —— **`Σ` を `n!` の素因数として取る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★原文の `Σ` の正体

原文は「ある有限素数集合 `Σ` の上で `ℤ[Σ⁻¹]` へ延びる」と書く。
★形式化では **`Σ = n!` の素因数集合**である——降下が段 `n` で立つとき、
`ℤ[1/n!]` が「`Σ` を反転した環」だからである。

★★本ファイルはその翻訳を与える:

    ch v ∉ (n!).primeFactors  ⟹  `n!` の冪は `v` の外

★★★これで `RatTowerLocalization.lean` の `conductorADiv_fin_eq_of_isUnit` の
仮定 `hM` が、原文どおりの「`Σ` の外の素点」から出る。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★★★`Σ` の外の素点では `n!` が可逆 -/

/-- ★★★★**`ch v` が `n!` の素因数でなければ `n!` は `v` の外**。

★機構は `hover`（`v` が `(ch v)` の上にある）——
`n! ∈ v` なら `ℤ` へ落として `ch v ∣ n!` となり、`ch v` が素数なので素因数になる。 -/
theorem fac_notMem_of_notMem_primeFactors (n : ℕ)
    (ch : FinitePlace F → ℕ)
    (hchprime : ∀ v, (ch v).Prime)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)}))
    (v : FinitePlace F) (hv : ch v ∉ (Nat.factorial n).primeFactors) :
    ((Nat.factorial n : ℕ) : NumberField.RingOfIntegers F) ∉ v.asIdeal := by
  intro hmem
  have hunder : Ideal.span {((ch v : ℕ) : ℤ)}
      = v.asIdeal.under ℤ := (hover v).over
  have hZ : ((Nat.factorial n : ℕ) : ℤ) ∈ v.asIdeal.under ℤ := by
    show ((Nat.factorial n : ℕ) : ℤ) ∈ v.asIdeal.comap (algebraMap ℤ _)
    simpa using hmem
  rw [← hunder, Ideal.mem_span_singleton] at hZ
  have hdvd : ch v ∣ Nat.factorial n := by exact_mod_cast hZ
  exact hv (Nat.mem_primeFactors.2 ⟨hchprime v, hdvd, Nat.factorial_ne_zero n⟩)

/-- ★★★★★★**`Σ` の外の素点では `n!` の冪すべてが `v` の外**。

★これが `conductorADiv_fin_eq_of_isUnit` の仮定 `hM` そのものである。 -/
theorem powers_fac_notMem (n : ℕ)
    (ch : FinitePlace F → ℕ)
    (hchprime : ∀ v, (ch v).Prime)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)}))
    (v : FinitePlace F) (hv : ch v ∉ (Nat.factorial n).primeFactors) :
    ∀ m ∈ Submonoid.powers ((Nat.factorial n : ℕ) : NumberField.RingOfIntegers F),
      m ∉ v.asIdeal := by
  rintro m ⟨k, rfl⟩ hmem
  exact fac_notMem_of_notMem_primeFactors F n ch hchprime hover v hv
    (v.isPrime.mem_of_pow_mem k hmem)


/-! ## ★★★★★★★★★★到達点 —— `Remark 1.5.1` の後半 -/

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 の後半** ——
`log-cond` の BD-class は対 `(X_ℚ, D_ℚ)` だけに依る。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★入力は `DivisorDescent.lean` の `exists_pair_descent` が出す
「有限段 `n` での同型＋可換正方形」と、
点の対応 `ePt` がその同型と両立すること（`hcompat`）である。

★★★**定数は `∑_{q ∣ n!} log q`** ——原文の `∑_{q ∈ Σ} log q` そのものであり、
**点にも定義体にも依らない**。

★★★★★鎖の全体:
`exists_pair_descent`（降下）→ `comap_eq_of_square`（幾何）
→ `conductorADiv_fin_eq_of_isUnit`（導手）→ `abs_logCond_sub_le_sum_log`（Σ 上の寄与）。 -/
theorem remark_1_5_1_bdeq
    {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    {n : ℕᵒᵖ}
    (A : Type) [CommRing A] [Algebra (NumberField.RingOfIntegers F) A]
    [IsLocalization
      (Submonoid.powers ((Nat.factorial n.unop : ℕ) : NumberField.RingOfIntegers F)) A]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ)
    (ch : FinitePlace F → ℕ) (hchprime : ∀ v, (ch v).Prime)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)}))
    (ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'))
    (hcompat : ∀ xF,
      Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ ePt xF
        = liftPointToBc F A hinv f xF ≫ φ ≫
          pullback.snd (overRatTowerDiagram.obj n).hom f')
    (hI : ∀ xF, pullbackIdeal F iZ.ker xF ≠ 0)
    (hI' : ∀ xF, pullbackIdeal F iZ'.ker (ePt xF) ≠ 0) :
    BDeq (fun xF => logCond F iZ.ker xF) (fun xF => logCond F iZ'.ker (ePt xF)) := by
  refine ⟨∑ q ∈ (Nat.factorial n.unop).primeFactors, Real.log q, fun xF => ?_⟩
  dsimp only
  refine abs_logCond_sub_le_sum_log F iZ.ker iZ'.ker xF (ePt xF) (hI xF) (hI' xF)
    (Nat.factorial n.unop).primeFactors
    (fun q hq => Nat.prime_of_mem_primeFactors hq) ch ?_ hover
  intro v hv
  exact conductorADiv_fin_eq_of_isUnit F _ A v
    (powers_fac_notMem F n.unop ch hchprime hover v hv) f f' iZ iZ' hinv φ ψ hsq
    xF (ePt xF) (hcompat xF) (hI xF) (hI' xF)

/-! ### ★出典の紐付け(`.src`) -/

def remark_1_5_1_bdeq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(後半——log-cond の BD-class が対 (X_ℚ, D_ℚ) だけに依ること)",
    sectionId := "genell-rem-1-5-1" }

def remark_1_5_1_bdeq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_pair_descent(対 (X, D) の spreading out)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pair_descent") 9,
    .citation "[ABC3]" "conductorADiv_fin_eq_of_isUnit(Σ の外で導手が一致)"
      (.inProject "ABC3" "ABC3.Found.GenEll.conductorADiv_fin_eq_of_isUnit") 9,
    .citation "[ABC3]" "abs_logCond_sub_le_sum_log(Σ 上の寄与の一様な上界)"
      (.inProject "ABC3" "ABC3.Found.GenEll.abs_logCond_sub_le_sum_log") 9,
    .implicitStep
      ("★★★★原文の `Σ` の正体は **`n!` の素因数集合**である" ++
       "——降下が段 n で立つとき、ℤ[1/n!] が「Σ を反転した環」だからである") 9,
    .implicitStep
      ("★★定数は ∑_{q ∣ n!} log q であり、**点にも定義体にも依らない**。" ++
       "原文が `≈`(BD-同値)で済ませているのはこの定数分である") 9,
    .implicitStep
      ("★★★残るもの: 点の対応 `ePt` とその両立 `hcompat` を**データとして受けている**。" ++
       "★原文では X_ℚ ≅ X′_ℚ が ℚ̄-点の対応を与え、" ++
       "𝒪_F へ戻すのに X′ の固有性(付値判定法)を使う。" ++
       "★★mathlib には固有射の付値判定法があるので、欠落ではなく我々の作業である") 9 ]

end ABC3.Found.GenEll
