import ABC3.Found.Falt1.AlmostDerivation

/-!
# [Falt1] Chapter I §4 —— `Theorem 4.1(ii)` の完全列(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §4、
`Theorem 4.1`(物理 p.13-14 = 印字 p.266-267)。

> *4.1. Theorem. (i) The map `Ω_{R/V} → Ω̄_{R̄/V}` induces almost
> isomorphisms `Ω_{R/V} ⊗_R R̄[1/p] ≅ Ω̄_{R̄/V}` and
> `Ω_{R/V} ⊗_R (R̄[1/p]/R̄) ≅ Ω̄_{R̄/R⊗_V V̄}`.
> (ii) The sequence
> `0 → Ω_{V̄/V} ⊗_{V̄} R̄ → Ω̄_{R̄/R} → Ω̄_{R̄/R⊗_V V̄} → 0` is exact.*
>
> *Proof. We only have to show that the first map of the exact sequence in
> (ii) is injective. We already know that its kernel is annihilated by `m`.
> But `Ω_{V̄/V} ⊗_V R̄` has no `m`-torsion.*

## 何が起きているか(2026-09-05 に判明したロードマップの訂正)

`falt1-goal.md` はこれまで「§3・§4 は §2 の壁(`Theorem 2.2`/`2.3`)を
継承する」と記録していたが、原典を逐語で読み直すと**そうではない**:

* `Theorem 3.1`/`3.2` の証明が使うのは `Definition 2.1` と
  **`Theorem 1.2`**(および局所コホモロジー・反射的加群)であって、
  `2.2`/`2.3` ではない。
* **`Theorem 4.1` の証明は上に引いた 3 行がすべて**であり、その入力は
  **`Theorem 2.4(i)` の核側**(「核は `m` で零化される」)である。

`Theorem 2.4(i)` は本トラックで証明済み(`AlmostDerivation.thm_2_4_i`)
なので、`Theorem 4.1(ii)` はそこから直ちに出る:

* **単射性**——核が `p^{2n}` で零化される(`thm_2_4_i` の核側)ことと、
  源 `Ω[A⁄R] ⊗_A B` に `p^{2n}`-捩れが無いことを合わせる
  (原文の *"But `Ω_{V̄/V} ⊗_V R̄` has no `m`-torsion"*)。
* **中央の完全性・右の全射性**——Jacobi–Zariski 完全列(mathlib の
  `KaehlerDifferential.exact_mapBaseChange_map` と `map_surjective`)。

## `Theorem 4.1(i)` について

★**訂正(2026-09-05)**: 以前ここに「(i) は `Ω̄`(`p` 進完備化した微分)を
含むので完備化の枠組みが要る」と書いていたが、**誤りだった**——原典を
260dpi で描画して逐語確認した結果、bar は `Ω` ではなく**添字の
`R`・`V`** に付いており(`Ω_{R̄/V̄}`)、`p` 進完備化ではない。
必要なのは通常の Kähler 微分と局所化 `R̄[1/p]` だけである。

その抽象的な内容 —— *"almost 同型は `[1/p]` を取ると本物の同型になる"*
—— は `thm_4_1_i_localized` で証明した(`localizedModule_map_bijective`
＋ `thm_2_4_i`)。残るのは `R∞`(単元の `p` 冪根を添加した塔)の設定と
`dlog(uᵢ)` を基底とする明示計算である。
-/

namespace ABC3.Found.Falt1

universe u

open KaehlerDifferential in
/-- **`Theorem 4.1(ii)`(Faltings, Ch.I §4)——完全列**。

`0 → Ω[A⁄R] ⊗_A B → Ω[B⁄R] → Ω[B⁄A] → 0`

原文の証明そのもの:**核が `m` で零化される**(`thm_2_4_i` の核側)
ことと**源に `m`-捩れが無い**ことを合わせて単射性が出る。中央の完全性と
右の全射性は Jacobi–Zariski から形式的に従う。 -/
theorem thm_4_1_ii {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    (p : A) (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p)))
    (hnotors : ∀ x : TensorProduct A B Ω[A⁄R], ((p^n) * (p^n)) • x = 0 → x = 0) :
    Function.Injective (KaehlerDifferential.mapBaseChange R A B) ∧
    Function.Exact (KaehlerDifferential.mapBaseChange R A B)
      (KaehlerDifferential.map R A B B) ∧
    Function.Surjective (KaehlerDifferential.map R A B B) := by
  refine ⟨?_, KaehlerDifferential.exact_mapBaseChange_map R A B,
    KaehlerDifferential.map_surjective R A B⟩
  rw [injective_iff_map_eq_zero]
  intro x hx
  exact hnotors x ((thm_2_4_i (R := R) p hAE hf0inj n w hw).1 x hx)

open KaehlerDifferential in
/-- **★`Theorem 4.1(i)` の抽象的な内容**——`Theorem 2.4(i)` の almost 同型は
**水準を反転すると本物の同型になる**。

原文の *"The map `Ω_{R/V} ⊗_R R̄ → Ω_{R̄/V̄}` induces almost isomorphisms
`Ω_{R/V} ⊗_R R̄[1/p] ≅ …`"* の機構そのもの——`[1/p]` を取ると
`m`-捩れが消え、almost 同型が honest な同型になる。

★**完備化の枠組みは要らない**(2026-09-05 に原典を 260dpi で再確認:
bar は `Ω` ではなく添字に付いており、`p` 進完備化ではない)。
局所化は完全なので**平坦性の議論も要らない**。 -/
theorem thm_4_1_i_localized {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    (p : A) (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p))) :
    Function.Bijective (LocalizedModule.map
      (Submonoid.powers (algebraMap A B (p ^ n * p ^ n)))
      (KaehlerDifferential.mapBaseChange R A B)) := by
  obtain ⟨hk, hc⟩ := thm_2_4_i (R := R) p hAE hf0inj n w hw
  refine localizedModule_map_bijective _ _ (fun ξ hξ => ?_) (fun x => ?_)
  · rw [algebraMap_smul]
    exact hk ξ hξ
  · rw [map_mul, mul_smul]
    exact Submodule.smul_mem _ _ (hc x)

/-! ## ★§4(b) の「`dlog(uᵢ)` を基底とする自由加群」(2026-09-05)

原文 §4(b):

> *The direct sums `R̄^d`, `R̄[1/p]^d`, etc. are free modules with the
> `dlog(uᵢ)` as basis.*

その根拠は「`V[T₁…T_d] → R` がエタール」という §4(a) の設定である。
一般に**形式的エタールなら微分は底変換で同型**になり、多項式環の微分は
自由(mathlib の `KaehlerDifferential.mvPolynomialBasis`)なので、
`Ω[R⁄V]` は `d` 個の `d(uᵢ)` を基底とする自由加群になる。 -/

/-- **形式的エタールなら微分は底変換で同型**——Jacobi–Zariski の
2 本の完全列(`H1Cotangent` 側と `Ω` 側)で両端が消えることから。 -/
theorem mapBaseChange_bijective_of_formallyEtale (R S T : Type u) [CommRing R] [CommRing S]
    [CommRing T] [Algebra R S] [Algebra S T] [Algebra R T] [IsScalarTower R S T]
    [Algebra.FormallyEtale S T] :
    Function.Bijective (KaehlerDifferential.mapBaseChange R S T) := by
  constructor
  · rw [injective_iff_map_eq_zero]
    intro x hx
    obtain ⟨y, hy⟩ := ((Algebra.H1Cotangent.exact_δ_mapBaseChange R S T) x).mp hx
    rw [← hy, Subsingleton.elim y 0, map_zero]
  · intro y
    exact ((KaehlerDifferential.exact_mapBaseChange_map R S T) y).mp (Subsingleton.elim _ _)

/-- **`Ω[R⁄V]` は `d(uᵢ)` を基底とする自由加群**——`V[Tᵢ] → R` が
形式的エタールなら。原文 §4(b) の *"free modules with the `dlog(uᵢ)`
as basis"* の内容(`uᵢ` が単元なら `d(uᵢ)` と `dlog(uᵢ)` は単元倍で
移り合う)。 -/
noncomputable def kaehlerBasisOfEtale (V : Type u) [CommRing V] (σ : Type u) (R : Type u)
    [CommRing R] [Algebra V R] [Algebra (MvPolynomial σ V) R]
    [IsScalarTower V (MvPolynomial σ V) R]
    [Algebra.FormallyEtale (MvPolynomial σ V) R] :
    Module.Basis σ R (Ω[R⁄V]) :=
  ((KaehlerDifferential.mvPolynomialBasis V σ).baseChange R).map
    (LinearEquiv.ofBijective (KaehlerDifferential.mapBaseChange V (MvPolynomial σ V) R)
      (mapBaseChange_bijective_of_formallyEtale V (MvPolynomial σ V) R))

/-- 基底の元はちょうど `d(uᵢ)`(`uᵢ := Tᵢ` の像)。 -/
theorem kaehlerBasisOfEtale_apply (V : Type u) [CommRing V] (σ : Type u) (R : Type u)
    [CommRing R] [Algebra V R] [Algebra (MvPolynomial σ V) R]
    [IsScalarTower V (MvPolynomial σ V) R]
    [Algebra.FormallyEtale (MvPolynomial σ V) R] (i : σ) :
    kaehlerBasisOfEtale V σ R i
      = KaehlerDifferential.D V R (algebraMap (MvPolynomial σ V) R (MvPolynomial.X i)) := by
  rw [kaehlerBasisOfEtale, Module.Basis.map_apply, Module.Basis.baseChange_apply,
    KaehlerDifferential.mvPolynomialBasis_apply]
  simp [KaehlerDifferential.mapBaseChange_tmul, KaehlerDifferential.map_D]

/-- 非空虚性——`V = ℤ`、`σ = Fin 1`、`R = Fin 2 → ℤ[X]`
(`ℤ[X]` 上分裂エタール、`Algebra.FormallyEtale.pi_iff`)。
`Ω[R⁄ℤ]` は階数 1 の自由 `R`-加群になる。 -/
noncomputable example :
    Module.Basis (Fin 1) (Fin 2 → MvPolynomial (Fin 1) ℤ)
      (Ω[(Fin 2 → MvPolynomial (Fin 1) ℤ)⁄ℤ]) := by
  haveI : Algebra.FormallyEtale (MvPolynomial (Fin 1) ℤ)
      (Fin 2 → MvPolynomial (Fin 1) ℤ) :=
    (Algebra.FormallyEtale.pi_iff (R := MvPolynomial (Fin 1) ℤ)
      (A := fun _ : Fin 2 => MvPolynomial (Fin 1) ℤ)).mpr (fun _ => inferInstance)
  exact kaehlerBasisOfEtale ℤ (Fin 1) (Fin 2 → MvPolynomial (Fin 1) ℤ)

/-! ## 非空虚性の対照

`R = ℤ`、`A = ℤ[X]`(**`Ω[A⁄R] ≠ 0`** なので主張が退化しない)、
`B = Fin 2 → A`(階数 2 の非自明な有限エタール拡大)、`p = 5`(真の非単元)。 -/

open KaehlerDifferential in
/-- 源 `B ⊗_A Ω[A⁄R] ≅ (Fin 2 → ℤ[X])` は `25`-捩れを持たない
(原文の *"`Ω_{V̄/V} ⊗_V R̄` has no `m`-torsion"* に当たる仮定)。 -/
theorem thm_4_1_nonvacuous_notors :
    ∀ x : TensorProduct (Polynomial ℤ) (Fin 2 → Polynomial ℤ) Ω[Polynomial ℤ⁄ℤ],
      (((5:Polynomial ℤ)^1) * ((5:Polynomial ℤ)^1)) • x = 0 → x = 0 := by
  intro x hx
  set e : TensorProduct (Polynomial ℤ) (Fin 2 → Polynomial ℤ) Ω[Polynomial ℤ⁄ℤ]
      ≃ₗ[Polynomial ℤ] (Fin 2 → Polynomial ℤ) :=
    (TensorProduct.congr (LinearEquiv.refl (Polynomial ℤ) (Fin 2 → Polynomial ℤ))
      (KaehlerDifferential.polynomialEquiv ℤ)).trans
      (TensorProduct.rid (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) with he
  have h1 : e ((((5:Polynomial ℤ)^1) * ((5:Polynomial ℤ)^1)) • x) = 0 := by rw [hx, map_zero]
  rw [map_smul] at h1
  have h2 : e x = 0 := by
    funext j
    have hj := congrFun h1 j
    simp only [Pi.smul_apply, Pi.zero_apply, smul_eq_mul] at hj
    rcases mul_eq_zero.mp hj with h | h
    · exact absurd h (by norm_num)
    · exact h
  exact e.injective (by rw [h2, map_zero])

theorem thm_4_1_nonvacuous_hf0inj :
    Function.Injective (algebraMap (Fin 2 → Polynomial ℤ)
      (Localization.Away (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ)
        (5 : Polynomial ℤ)))) := by
  refine IsLocalization.injective
    (M := Submonoid.powers
      (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ) (5 : Polynomial ℤ))) _ ?_
  rw [Submonoid.powers_le, mem_nonZeroDivisors_iff]
  have key : ∀ z : Fin 2 → Polynomial ℤ,
      (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) 5 * z = 0 → z = 0 := by
    intro z hz
    funext j
    have hj := congrFun hz j
    simp only [Pi.mul_apply, Pi.zero_apply] at hj
    rw [show (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) 5 j = 5 by simp] at hj
    rcases mul_eq_zero.mp hj with h | h
    · exact absurd h (by norm_num)
    · exact h
  exact ⟨key, fun z hz => key z (by rw [mul_comm]; exact hz)⟩

open KaehlerDifferential in
/-- **`thm_4_1_ii` の非空虚性**——仮定を全て満たす具体例が実在する。 -/
example :
    Function.Injective
      (KaehlerDifferential.mapBaseChange ℤ (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) ∧
    Function.Exact
      (KaehlerDifferential.mapBaseChange ℤ (Polynomial ℤ) (Fin 2 → Polynomial ℤ))
      (KaehlerDifferential.map ℤ (Polynomial ℤ) (Fin 2 → Polynomial ℤ)
        (Fin 2 → Polynomial ℤ)) ∧
    Function.Surjective (KaehlerDifferential.map ℤ (Polynomial ℤ) (Fin 2 → Polynomial ℤ)
      (Fin 2 → Polynomial ℤ)) := by
  have hAE : IsAlmostEtaleCovering (A := Polynomial ℤ) (B := Fin 2 → Polynomial ℤ)
      (5 : Polynomial ℤ) := isAlmostEtaleCovering_of_etale_general _
  refine thm_4_1_ii (R := ℤ) (5 : Polynomial ℤ) hAE thm_4_1_nonvacuous_hf0inj 1
    ((5 : Polynomial ℤ)^1 • Algebra.FormallyUnramified.elem
      (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) ?_ thm_4_1_nonvacuous_notors
  letI := awayAlgebra (5 : Polynomial ℤ) (A := Polynomial ℤ) (B := Fin 2 → Polynomial ℤ)
  haveI := hAE.2.2.1
  haveI := (hAE.2.1 : Module.Finite _ _)
  rw [map_smul, diagonalCompare_elem_eq]

open KaehlerDifferential in
/-- **`thm_4_1_i_localized` の非空虚性**——同じ具体例
(`R = ℤ`・`A = ℤ[X]`・`B = Fin 2 → A`・`p = 5`)で、`25` を反転すると
`Ω[A⁄R] ⊗_A B → Ω[B⁄R]` が**本物の全単射**になる。 -/
example : Function.Bijective (LocalizedModule.map
    (Submonoid.powers (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ)
      ((5 : Polynomial ℤ) ^ 1 * (5 : Polynomial ℤ) ^ 1)))
    (KaehlerDifferential.mapBaseChange ℤ (Polynomial ℤ) (Fin 2 → Polynomial ℤ))) := by
  have hAE : IsAlmostEtaleCovering (A := Polynomial ℤ) (B := Fin 2 → Polynomial ℤ)
      (5 : Polynomial ℤ) := isAlmostEtaleCovering_of_etale_general _
  refine thm_4_1_i_localized (R := ℤ) (5 : Polynomial ℤ) hAE thm_4_1_nonvacuous_hf0inj 1
    ((5 : Polynomial ℤ) ^ 1 • Algebra.FormallyUnramified.elem
      (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) ?_
  letI := awayAlgebra (5 : Polynomial ℤ) (A := Polynomial ℤ) (B := Fin 2 → Polynomial ℤ)
  haveI := hAE.2.2.1
  haveI := (hAE.2.1 : Module.Finite _ _)
  rw [map_smul, diagonalCompare_elem_eq]

/-! ## 原文どおりの形——「核は `m` で零化される。しかし源に `m`-捩れは無い」

`p`-可除な塔 `PDivTower` の上で述べると、原文の証明の 3 行

> *We already know that its kernel is annihilated by `m`.
> But `Ω_{V̄/V} ⊗_V R̄` has **no `m`-torsion**.*

がそのまま仮定と結論になる。 -/

open KaehlerDifferential in
/-- **`Theorem 4.1(ii)`、塔の上での形**。仮定 `hnotors` が原文の
*"has no `m`-torsion"* に、`thm_2_4_i_tower` が *"its kernel is
annihilated by `m`"* に対応する。 -/
theorem thm_4_1_ii_tower {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (hq : 2 ≤ q) (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0)))))
    (hnotors : ∀ (k : ℕ) (ξ : TensorProduct A B Ω[A⁄R]), T.ϖ k • ξ = 0 → ξ = 0) :
    Function.Injective (KaehlerDifferential.mapBaseChange R A B) ∧
    Function.Exact (KaehlerDifferential.mapBaseChange R A B)
      (KaehlerDifferential.map R A B B) ∧
    Function.Surjective (KaehlerDifferential.map R A B B) := by
  refine ⟨?_, KaehlerDifferential.exact_mapBaseChange_map R A B,
    KaehlerDifferential.map_surjective R A B⟩
  rw [injective_iff_map_eq_zero]
  intro x hx
  exact hnotors 0 x ((thm_2_4_i_tower (R := R) hq T hAET hf0inj).1 0 x hx)

end ABC3.Found.Falt1
