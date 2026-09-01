/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicEvalGen
import ABC3.Found.GaloisRep.TateSpecialize
import ABC3.Meta.Claim

/-!
# 第 1126 ブロック —— **係数環つきの評価も環準同型である**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か

第 96（`TateSpecialize.lean`）は **`ℤ⟦q⟧ →+* R`**（`evalAdicHom`）を与える。
第 1106（`AdicEvalGen.lean`）は係数環つきの評価 `evalAdicMap φ` を**関数として**与える。
★本ブロックはその 2 つを繋ぎ、**`S⟦q⟧ →+* R`**（`evalAdicMapHom φ q hq`）を立てる。

☆機構は第 96 とまったく同じ——打ち切り（`PowerSeries.trunc`）で部分和を多項式の
評価に書き換え、`trunc_trunc_mul_trunc` で法 `I^n` の乗法性を出す。
★`Int.castRingHom R` を `φ : S →+* R` に差し替えるだけである。

## ★★★★これで何が繋がるか

`Lemma 3.5` から `d + 1 < l` を外す道（§3 の枠の節点 3）の降下は

    `A₀ = PowerSeries ℤ[ζ_l]`  →（`PowerSeries.map`）→  `A₁ = PowerSeries ℚ(ζ_l)`

で行い、最後に `A₀ → R`（`ζ_l ↦ ζ`、`X ↦ q`）で特殊化する。
★その最後の段が本ブロックの `evalAdicMapHom` である。
☆第 1125（`map_*DF`）が「DF 形は環準同型を通り抜ける」を与えているので、
`evalAdicMapHom` を通せば `R` の中の等式になる。
-/

namespace ABC3.Found.GaloisRep

variable {R S : Type} [CommRing R] [CommRing S] {I : Ideal R}

/-! ## ★打ち切り -/

/-- ★係数環つきの部分和も「打ち切り多項式の評価」である。 -/
theorem partialEvalMap_eq_eval_trunc (φ : S →+* R) (f : PowerSeries S) (q : R) (n : ℕ) :
    partialEvalMap φ f q n = (PowerSeries.trunc n f).eval₂ φ q := by
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    have h0 : PowerSeries.trunc 0 f = 0 := by
      ext m
      rw [PowerSeries.coeff_trunc]
      simp
    rw [h0, partialEvalMap]
    simp
  · obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    rw [Polynomial.eval₂_eq_sum_range' φ (PowerSeries.natDegree_trunc_lt f m) q,
      partialEvalMap]
    refine Finset.sum_congr rfl (fun k hk => ?_)
    rw [PowerSeries.coeff_trunc, if_pos (Finset.mem_range.1 hk)]

/-- ★低次の係数が消える多項式は `q ∈ I` で `I^n` に値を取る（係数環つき）。 -/
theorem eval_mem_pow_of_coeff_zero_map (φ : S →+* R) (p : Polynomial S) (n : ℕ)
    (hp : ∀ k, k < n → p.coeff k = 0) {q : R} (hq : q ∈ I) : p.eval₂ φ q ∈ I ^ n := by
  rw [Polynomial.eval₂_eq_sum_range' φ (Nat.lt_succ_self p.natDegree) q]
  refine Submodule.sum_mem _ (fun k _ => ?_)
  rcases lt_or_ge k n with hk | hk
  · rw [hp k hk]
    simp
  · have hpow : q ^ k = q ^ n * q ^ (k - n) := by
      rw [← pow_add]
      congr 1
      omega
    rw [hpow, ← mul_assoc]
    exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n))

/-! ## ★★★乗法性 -/

/-- ★★★係数環つきの部分和も法 `I^n` で乗法的である。 -/
theorem partialEvalMap_mul_smodEq (φ : S →+* R) (f g : PowerSeries S) {q : R} (hq : q ∈ I)
    (n : ℕ) :
    partialEvalMap φ (f * g) q n ≡ partialEvalMap φ f q n * partialEvalMap φ g q n
      [SMOD (I ^ n • ⊤ : Submodule R R)] := by
  rw [SModEq.sub_mem, partialEvalMap_eq_eval_trunc, partialEvalMap_eq_eval_trunc,
    partialEvalMap_eq_eval_trunc, ← Polynomial.eval₂_mul, ← Polynomial.eval₂_sub]
  have hcoeff : ∀ k, k < n →
      (PowerSeries.trunc n (f * g)
        - (PowerSeries.trunc n f) * (PowerSeries.trunc n g)).coeff k = 0 := by
    intro k hk
    rw [Polynomial.coeff_sub, PowerSeries.coeff_trunc, if_pos hk]
    have h1 := PowerSeries.trunc_trunc_mul_trunc (n := n) f g
    have h2 := congrArg (fun p : Polynomial S => p.coeff k) h1
    simp only [PowerSeries.coeff_trunc, if_pos hk] at h2
    have h3 : Polynomial.coeff ((PowerSeries.trunc n f) * (PowerSeries.trunc n g)) k
        = PowerSeries.coeff k
            ((((PowerSeries.trunc n f) : Polynomial S) : PowerSeries S)
              * (((PowerSeries.trunc n g) : Polynomial S) : PowerSeries S)) := by
      rw [← Polynomial.coe_mul, Polynomial.coeff_coe]
    rw [h3, h2]
    ring
  have hmem := eval_mem_pow_of_coeff_zero_map φ
    (PowerSeries.trunc n (f * g) - (PowerSeries.trunc n f) * (PowerSeries.trunc n g)) n hcoeff hq
  simpa using hmem

/-! ## ★★★★環準同型の各条件 -/

/-- ★★★**係数環つきの評価は加法的である**。 -/
theorem evalAdicMap_add [IsAdicComplete I R] (φ : S →+* R) (f g : PowerSeries S)
    (q : R) (hq : q ∈ I) :
    evalAdicMap φ (f + g) q hq = evalAdicMap φ f q hq + evalAdicMap φ g q hq := by
  refine evalAdicMap_unique φ _ _ _ _ (fun n => ?_)
  have hps : partialEvalMap φ (f + g) q n
      = partialEvalMap φ f q n + partialEvalMap φ g q n := by
    simp only [partialEvalMap, map_add, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl (fun k _ => by ring)
  rw [hps]
  exact SModEq.add (evalAdicMap_spec φ f q hq n) (evalAdicMap_spec φ g q hq n)

/-- ★★★**係数環つきの評価は乗法的である**。 -/
theorem evalAdicMap_mul [IsAdicComplete I R] (φ : S →+* R) (f g : PowerSeries S)
    (q : R) (hq : q ∈ I) :
    evalAdicMap φ (f * g) q hq = evalAdicMap φ f q hq * evalAdicMap φ g q hq := by
  refine evalAdicMap_unique φ _ _ _ _ (fun n => ?_)
  refine (partialEvalMap_mul_smodEq φ f g hq n).trans ?_
  have ha := evalAdicMap_spec φ f q hq n
  have hb := evalAdicMap_spec φ g q hq n
  rw [SModEq.sub_mem] at ha hb ⊢
  have hexp : partialEvalMap φ f q n * partialEvalMap φ g q n
      - evalAdicMap φ f q hq * evalAdicMap φ g q hq
      = (partialEvalMap φ f q n - evalAdicMap φ f q hq) * partialEvalMap φ g q n
        + evalAdicMap φ f q hq * (partialEvalMap φ g q n - evalAdicMap φ g q hq) := by ring
  rw [hexp]
  refine Submodule.add_mem _ ?_ ?_
  · simpa using Ideal.mul_mem_right _ _ (by simpa using ha)
  · simpa using Ideal.mul_mem_left _ _ (by simpa using hb)

/-- ★定数の値。 -/
theorem evalAdicMap_C [IsAdicComplete I R] (φ : S →+* R) (c : S) (q : R) (hq : q ∈ I) :
    evalAdicMap φ (PowerSeries.C c) q hq = φ c := by
  refine evalAdicMap_unique φ _ _ _ _ (fun n => ?_)
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    have htop : (I ^ 0 • (⊤ : Submodule R R)) = ⊤ := by simp
    rw [SModEq.sub_mem, htop]
    exact Submodule.mem_top
  · have hp : partialEvalMap φ (PowerSeries.C c) q n = φ c := by
      obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
      rw [partialEvalMap, Finset.sum_range_succ']
      have hz : ∀ x ∈ Finset.range m,
          φ (PowerSeries.coeff (x + 1) (PowerSeries.C c)) * q ^ (x + 1) = 0 := by
        intro x _
        rw [PowerSeries.coeff_C, if_neg (Nat.succ_ne_zero x)]
        simp
      rw [Finset.sum_eq_zero hz]
      simp
    rw [hp]

theorem evalAdicMap_one [IsAdicComplete I R] (φ : S →+* R) (q : R) (hq : q ∈ I) :
    evalAdicMap φ (1 : PowerSeries S) q hq = 1 := by
  have h1 : (1 : PowerSeries S) = PowerSeries.C 1 := by simp
  rw [h1, evalAdicMap_C]
  simp

/-- ★★★★★★★★**係数環つきの評価は環準同型である** `S⟦q⟧ →+* R`（第 1126）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`φ : S →+* R` と `q ∈ I` から `S⟦q⟧ → R`（`X ↦ q`、係数は `φ`）。
☆`S = ℤ`・`φ = Int.castRingHom R` なら第 96 の `evalAdicHom` に一致する。 -/
noncomputable def evalAdicMapHom [IsAdicComplete I R] (φ : S →+* R) (q : R) (hq : q ∈ I) :
    PowerSeries S →+* R where
  toFun := fun f => evalAdicMap φ f q hq
  map_one' := evalAdicMap_one φ q hq
  map_mul' := fun f g => evalAdicMap_mul φ f g q hq
  map_zero' := by
    have h0 : (0 : PowerSeries S) = PowerSeries.C 0 := by simp
    have hc := evalAdicMap_C (I := I) φ 0 q hq
    rw [← h0] at hc
    simpa using hc
  map_add' := fun f g => evalAdicMap_add φ f g q hq

@[simp] theorem evalAdicMapHom_apply [IsAdicComplete I R] (φ : S →+* R) (q : R) (hq : q ∈ I)
    (f : PowerSeries S) : evalAdicMapHom φ q hq f = evalAdicMap φ f q hq := rfl

/-- ★★★★**`S = ℤ` なら第 96 の `evalAdicHom` に一致する**。 -/
theorem evalAdicMapHom_int [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    evalAdicMapHom (Int.castRingHom R) q hq = evalAdicHom q hq := by
  ext f
  exact evalAdicMap_int f q hq

/-! ## ★出典の紐付け(`.src`) -/

def evalAdicMapHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(係数環つきの q-展開の評価が環準同型であること)",
    sectionId := "genell-def-3-3" }

def evalAdicMapHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "PowerSeries.trunc_trunc_mul_trunc(打ち切りは n 次未満で乗法的)"
      (.inMathlib "PowerSeries.trunc_trunc_mul_trunc") 1,
    .citation "[ABC3]" "evalAdicHom(係数が ℤ の場合、第 96、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdicHom") 1,
    .implicitStep
      ("★★**2026-09-01（第 1126）**——第 96 の証明の `Int.castRingHom R` を " ++
       "`φ : S →+* R` に差し替えただけである。" ++
       "☆これで `PowerSeries ℤ[ζ_l] →+* R`（`ζ_l ↦ ζ`、`X ↦ q`）が立ち、" ++
       "第 1125 の `map_*DF` と合わせて節点 3 の降下路が繋がる。") 3 ]

end ABC3.Found.GaloisRep
