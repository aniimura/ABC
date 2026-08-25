/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm64Slim
import ABC3.Found.Divisor.Ex63Datum
import ABC3.Found.NumberField.GaloisFaithfulAll
import ABC3.Found.Divisor.ArithRlfFlat

/-!
# [FrdI] Theorem 6.4, (i) —— 算術の `𝒟` は `Φ` に関して Div-slim（圏側の配線）

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.115。

原文 (FrdI p.115):
> that D is Div-slim [relative to Φ, hence also relative to Φpf, Φrlf].

## ★★何をする節点か

`Thm64Slim.lean` は **`Φ` から `Φ^pf` / `Φ^rlf` への伝播**（角括弧の `hence also`）を作った。
`GaloisFaithfulBase.lean` は **数論の側**（数体の自己同型で素イデアルを固定するものは恒等、
★Chebotarev を使わない）を作った。★本ファイルはその 2 つを**繋ぐ圏側の配線**である。

## ★★★3 段の内訳

| 段 | 主張 | 中身 |
|---|---|---|
| 1 | `isDivSlim_of_endo_map_injective` | `Aut(𝒟_A → 𝒟)` の成分は `X.left` の**自己準同型**なので、`Φ` が自己準同型を区別すれば Div-slim |
| 2 | `resPlaceOf_eq_of_pullOf_eq` | 引き戻しが有効因子の上で一致 ⟹ **素点の制限が一致**（`localDeg = 1` ＋ 各素点に台をもつ正の因子） |
| 3 | `finSubOp_isDivSlim` | 素点の制限が一致 ⟹ 射が一致（★数論の側） |

★★段 2 の要点は **`FinSub` の自己射は全単射**（有限次拡大の体の自己準同型）なので
`arithExtend` の局所次数因子が `1` になり（`localDeg_eq_one_of_bijective'`）、
引き戻しが「素点を引き戻すだけ」になることである。

★そこから各素点 `v` に台をもつ**正の**算術因子（`isGenSubgroupR_arithDivGroup`）を
当てれば、`resPlace` の一致が読み取れる。

## ★★★★訂正の記録 —— [Mzk7] Cor 1.1.6 は要らない

当初この節点には [Mzk7] `Corollary 1.1.6` が要ると見積もっていたが、
★`Φ.map e` が**素点の引き戻し**なので、Div-slim は
「`resPlace` が射を区別する」＝**数論の言明**に帰着し、圏論の入力は要らない。
（`Frobenius-slim` の側は依然 [Mzk7] `Corollary 1.1.6` を使う —— そちらは
`Aut (Over.forget A)` を `Gal(K̄/F)` に埋める必要があるためである。）
-/

namespace ABC3.Found.FrdI

open CategoryTheory ABC3.Found.Divisor _root_.NumberField
open scoped _root_.NumberField NNReal

universe v u w

/-! ## ★段 1. 一般形 —— 自己準同型を区別すれば Div-slim -/

/-- ★★★★**`Φ` が自己準同型を区別すれば `𝒟` は `Φ` に関して Div-slim**。

★`overPhiAut Φ A η` の成分は `Φ.functor.map (η.hom.app X).op` で、
`η.hom.app X : X.left ⟶ X.left` は **`𝒟` の自己準同型**である。
したがって `Φ` が自己準同型の上で単射なら、`overPhiAut Φ A` は単射になる。

★★**`𝒞` を一切見ない**のが `Definition 4.5, (iv)` の特徴で、ここでもそれが効く。 -/
theorem isDivSlim_of_endo_map_injective {D : Type u} [Category.{v} D]
    {Φ : MonoidOn.{v, u, w} D}
    (hfaith : ∀ (X : D) (f g : X ⟶ X), Φ.map f = Φ.map g → f = g) :
    IsDivSlim Φ := by
  intro A η₁ η₂ hη
  apply Iso.ext
  apply NatTrans.ext
  funext X
  have hX : (overPhiAut Φ A η₁).hom.app (Opposite.op X)
      = (overPhiAut Φ A η₂).hom.app (Opposite.op X) := by rw [hη]
  exact hfaith X.left _ _ (congrArg AddCommMonCat.Hom.hom hX)

/-- ★★★**`ArithDatum` の水準の Div-slim** —— 引き戻しが
**有効な**算術因子の上で自己準同型を区別すればよい。

★`Φ(A) = Γ_A ∩ ℝ≥0[V_A]` なので、`Φ.map` の一致は
「群 `Γ_A` の**有効**な元の上での引き戻しの一致」に他ならない。 -/
theorem arithDatum_isDivSlim {D : Type u} [Category.{v} D]
    (Δ : ArithDatum.{v, u, w} D) (hD : IsOfFSMType D)
    (hfaith : ∀ (X : D) (f g : X ⟶ X),
      (∀ x : Δ.primes X →₀ ℝ, x ∈ Δ.grp X → (∀ s, 0 ≤ x s) → Δ.pull f x = Δ.pull g x) →
      f = g) :
    IsDivSlim (Δ.phi hD) := by
  refine isDivSlim_of_endo_map_injective (fun X f g hfg => hfaith X f g ?_)
  intro x hx hxpos
  have hmem : x ∈ effR (Δ.grp X) := mem_effR.mpr ⟨hx, hxpos⟩
  exact congrArg (fun t : effR (Δ.grp X) →+ effR (Δ.grp X) =>
    (t ⟨x, hmem⟩ : Δ.primes X →₀ ℝ)) hfg

/-! ## ★段 2. 引き戻しの一致から素点の制限の一致へ -/

section Arith

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]

/-- ★**射 `f : L ⟶ M` に沿った素点の制限** `V(M) → V(L)`。 -/
noncomputable def resPlaceOf {L M : FinSub F Kbar} (f : L ⟶ M) :
    ArithPlace M.toIF → ArithPlace L.toIF :=
  letI := algOfHom f
  resPlace

/-- ★★**全単射な射に沿った引き戻しは「素点を引き戻すだけ」**。

★局所次数の因子が `1` になる（`localDeg_eq_one_of_bijective'`）。 -/
theorem pullOf_apply_of_bijective {L M : FinSub F Kbar} (f : L ⟶ M)
    (hb : Function.Bijective (FinSub.hom f))
    (x : ArithPlace L.toIF →₀ ℝ) (V : ArithPlace M.toIF) :
    pullOf f x V = x (resPlaceOf f V) := by
  letI := algOfHom f
  exact arithExtend_apply' hb x V

/-- ★★★★**有効な算術因子の引き戻しが一致すれば素点の制限が一致する**。

★★中身は「各素点 `v` に台をもつ**正の**算術因子がある」
（`isGenSubgroupR_arithDivGroup`）だけである ——
`v := resPlaceOf f V` に台をもつ因子を当てれば、
`resPlaceOf g V ≠ v` なら右辺が `0` になって矛盾する。

★`FinSub` の自己射が全単射であること（有限次拡大の体の自己準同型）を使う。 -/
theorem resPlaceOf_eq_of_pullOf_eq {L : FinSub F Kbar} (f g : L ⟶ L)
    (h : ∀ x : ArithPlace L.toIF →₀ ℝ, x ∈ arithDivGroup L.toIF → (∀ s, 0 ≤ x s) →
      pullOf f x = pullOf g x) :
    resPlaceOf f = resPlaceOf g := by
  classical
  haveI := L.fin
  have hbf : Function.Bijective (FinSub.hom f) := (FinSub.hom f).bijective
  have hbg : Function.Bijective (FinSub.hom g) := (FinSub.hom g).bijective
  funext V
  set v : ArithPlace L.toIF := resPlaceOf f V with hv
  obtain ⟨r, hr, hmem⟩ := (isGenSubgroupR_arithDivGroup (L := L.toIF)) v
  have hnn : ∀ s, 0 ≤ (Finsupp.single v r) s := by
    intro s
    rcases eq_or_ne v s with rfl | hvs
    · simpa using hr.le
    · simp [hvs]
  have h2 := congrArg (fun t : ArithPlace L.toIF →₀ ℝ => t V)
    (h (Finsupp.single v r) hmem hnn)
  simp only [pullOf_apply_of_bijective f hbf, pullOf_apply_of_bijective g hbg] at h2
  rw [← hv] at h2
  simp only [Finsupp.single_eq_same] at h2
  by_contra hne
  rw [Finsupp.single_apply, if_neg hne] at h2
  exact hr.ne' h2

/-! ## ★段 3. 算術の `𝒟` は Div-slim -/

variable (F Kbar) in
/-- ★★★★★**[FrdI] Theorem 6.4, (i)** —— 算術の `𝒟 = B(G)⁰` は
`Φ` に関して **Div-slim**。

原文 (FrdI p.115):
> that D is Div-slim [relative to Φ, hence also relative to Φpf, Φrlf].

★★仮定 `hfaith` が**数論の言明**（素点の制限が一致する自己射は等しい）で、
その中身は `GaloisFaithfulBase.lean` の
`eq_one_of_fixes_prime_of_galoisClosure`（★Chebotarev を使わない）である。

★★★角括弧の `hence also`（`Φ^pf` / `Φ^rlf` への伝播）は
`isDivSlim_pfOn` / `isDivSlim_phiScOn_nnreal`（`Thm64Slim.lean`）が与える。 -/
theorem finSubOp_isDivSlim [IsGalois F Kbar]
    (hfaith : ∀ (L : FinSub F Kbar) (f g : L ⟶ L), resPlaceOf f = resPlaceOf g → f = g) :
    IsDivSlim (phiGalois F Kbar) := by
  refine arithDatum_isDivSlim (arithDatumGalois F Kbar) finSubOp_isOfFSMType ?_
  intro X f g hfg
  refine Quiver.Hom.unop_inj (hfaith X.unop f.unop g.unop ?_)
  exact resPlaceOf_eq_of_pullOf_eq f.unop g.unop (fun x hx hxpos => hfg x hx hxpos)


/-! ## ★段 3 の中身 —— 素点の制限から射へ（数論の側） -/

/-- ★**`𝓞` の水準の射**（★一般の `[Algebra L M]` の側で固定する）。

★★`L = M` の文脈で `algebraMap (𝓞 L) (𝓞 L)` と書き直すと
`Algebra.id` を拾って**恒等写像**になってしまう（`tools\lean-idioms.md`）。
ここで一度だけ固定しておく。 -/
noncomputable def intMapOf {L M : FinSub F Kbar} (f : L ⟶ M) : 𝓞 L.toIF →+* 𝓞 M.toIF :=
  letI := algOfHom f
  algebraMap (𝓞 L.toIF) (𝓞 M.toIF)

/-- ★`intMapOf` は体の側では `f` そのもの。 -/
theorem coe_intMapOf {L M : FinSub F Kbar} (f : L ⟶ M) (x : 𝓞 L.toIF) :
    ((intMapOf f x : 𝓞 M.toIF) : M.toIF) = (FinSub.hom f) (x : L.toIF) := rfl

/-- ★★**全単射な射は `𝓞` の水準でも全単射**。 -/
theorem intMapOf_bijective {L M : FinSub F Kbar} (f : L ⟶ M)
    (hb : Function.Bijective (FinSub.hom f)) : Function.Bijective (intMapOf f) := by
  letI := algOfHom f
  refine ⟨fun x y hxy => RingOfIntegers.coe_injective
    (hb.1 (by rw [← coe_intMapOf f x, ← coe_intMapOf f y, hxy])), ?_⟩
  exact ringOfIntegers_surjective_of_bijective (L := L.toIF) (M := M.toIF) hb

/-- ★★★**素点の制限の一致は素イデアルの引き戻しの一致**。

★有限素点の成分を `FinitePlace.mk V` で読むだけである
（`FinitePlace.maximalIdeal_mk` が `HeightOneSpectrum` を復元する）。 -/
theorem comap_intMapOf_eq {L : FinSub F Kbar} (f g : L ⟶ L) (h : resPlaceOf f = resPlaceOf g)
    (V : IsDedekindDomain.HeightOneSpectrum (𝓞 L.toIF)) :
    V.asIdeal.comap (intMapOf f) = V.asIdeal.comap (intMapOf g) := by
  have h1 := congrFun h (Sum.inl (FinitePlace.mk V))
  simp only [resPlaceOf, resPlace_inl, resFin, FinitePlace.maximalIdeal_mk, Sum.inl.injEq] at h1
  have h2 := congrArg NumberField.FinitePlace.maximalIdeal h1
  simp only [FinitePlace.maximalIdeal_mk] at h2
  exact congrArg IsDedekindDomain.HeightOneSpectrum.asIdeal h2

open scoped Pointwise in
/-- ★★★★★★**素点の制限が一致する自己射は等しい** —— `Theorem 6.4, (i)` の
Div-slim の**数論の側**（★仮定なし）。

原文 (FrdI p.115):
> since any automorphism of a number field that fixes all of the valuations of the
> number field is clearly equal to the identity automorphism

★★★中身:

| 段 | 根拠 |
|---|---|
| `resPlace` の一致 → 素イデアルの引き戻しの一致 | `comap_intMapOf_eq` |
| `τ := g⁻¹ ∘ f` は全素イデアルを固定する | `intMapOf g` が全単射（引き戻しの像を動かす） |
| `τ = 1` | ★`eq_one_of_fixes_all_primes`（Chebotarev 不使用） |
| `f = g` | `g (τ x) = f x` に `τ = 1` を入れる | -/
theorem finSub_endo_eq_of_resPlaceOf_eq {L : FinSub F Kbar} (f g : L ⟶ L)
    (h : resPlaceOf f = resPlaceOf g) : f = g := by
  classical
  haveI := L.fin
  have hbf : Function.Bijective (FinSub.hom f) := (FinSub.hom f).bijective
  have hbg : Function.Bijective (FinSub.hom g) := (FinSub.hom g).bijective
  set eg : L.toIF ≃ₐ[F] L.toIF := AlgEquiv.ofBijective (FinSub.hom g) hbg with heg
  set ef : L.toIF ≃ₐ[F] L.toIF := AlgEquiv.ofBijective (FinSub.hom f) hbf with hef
  set ρ : L.toIF ≃ₐ[F] L.toIF := ef.trans eg.symm with hρ
  set τ : L.toIF ≃ₐ[ℚ] L.toIF := { ρ.toRingEquiv with commutes' := fun q => by simp } with hτdef
  have hgτ : ∀ x : L.toIF, (FinSub.hom g) (τ x) = (FinSub.hom f) x := by
    intro x
    show (FinSub.hom g) (eg.symm ((FinSub.hom f) x)) = _
    exact eg.apply_symm_apply _
  have hbφg : Function.Bijective (intMapOf g) := intMapOf_bijective g hbg
  set τO : 𝓞 L.toIF →+* 𝓞 L.toIF :=
    MulSemiringAction.toRingHom (L.toIF ≃ₐ[ℚ] L.toIF) (𝓞 L.toIF) τ with hτO
  have hcomp : (intMapOf g).comp τO = intMapOf f := by
    refine RingHom.ext (fun x => RingOfIntegers.coe_injective ?_)
    show (FinSub.hom g) (τ (algebraMap (𝓞 L.toIF) L.toIF x))
      = (FinSub.hom f) (algebraMap (𝓞 L.toIF) L.toIF x)
    exact hgτ _
  have hsurj : Function.Surjective τO := by
    intro y
    refine ⟨τ.symm • y, RingOfIntegers.coe_injective ?_⟩
    show τ (τ.symm (algebraMap (𝓞 L.toIF) L.toIF y)) = algebraMap (𝓞 L.toIF) L.toIF y
    exact τ.apply_symm_apply _
  have hτ1 : τ = 1 := by
    refine ABC3.Found.NF.eq_one_of_fixes_all_primes L.toIF τ ?_
    intro 𝔭 hp hne
    haveI := hp
    set eφg : 𝓞 L.toIF ≃+* 𝓞 L.toIF := RingEquiv.ofBijective (intMapOf g) hbφg with heφg
    set q : Ideal (𝓞 L.toIF) := 𝔭.comap (eφg.symm : 𝓞 L.toIF →+* 𝓞 L.toIF) with hq
    haveI hqp : q.IsPrime := Ideal.comap_isPrime _ _
    have hgq : q.comap (intMapOf g) = 𝔭 := by
      rw [hq, Ideal.comap_comap]
      have hid : (eφg.symm : 𝓞 L.toIF →+* 𝓞 L.toIF).comp (intMapOf g) = RingHom.id _ :=
        RingHom.ext (fun y => eφg.symm_apply_apply y)
      rw [hid, Ideal.comap_id]
    have hqne : q ≠ ⊥ := by
      intro h0
      rw [h0, Ideal.comap_bot_of_injective _ hbφg.1] at hgq
      exact hne hgq.symm
    have hV := comap_intMapOf_eq f g h ⟨q, hqp, hqne⟩
    have hfq : q.comap (intMapOf f) = 𝔭.comap τO := by
      rw [← hcomp, ← Ideal.comap_comap, hgq]
    rw [hfq, hgq] at hV
    have hmc := Ideal.map_comap_of_surjective τO hsurj 𝔭
    rw [hV] at hmc
    rw [Ideal.pointwise_smul_def]
    exact hmc
  refine AlgHom.ext (fun x => ?_)
  have hx := congrArg (fun t : L.toIF ≃ₐ[ℚ] L.toIF => t x) hτ1
  simp only [AlgEquiv.one_apply] at hx
  have h2 : (FinSub.hom g) (τ x) = (FinSub.hom f) x := hgτ x
  rw [hx] at h2
  exact h2.symm

variable (F Kbar) in
/-- ★★★★★★**[FrdI] Theorem 6.4, (i)** —— 算術の `𝒟 = B(G)⁰` は
`Φ` に関して **Div-slim**（★仮定なし）。

原文 (FrdI p.115):
> that D is Div-slim [relative to Φ, hence also relative to Φpf, Φrlf].

★★★角括弧の `hence also` は `isDivSlim_pfOn` / `isDivSlim_phiScOn_nnreal`
（`Thm64Slim.lean`）が与える。 -/
theorem finSubOp_isDivSlim' [IsGalois F Kbar] : IsDivSlim (phiGalois F Kbar) :=
  finSubOp_isDivSlim F Kbar (fun _ f g hfg => finSub_endo_eq_of_resPlaceOf_eq f g hfg)


/-! ## ★★★★角括弧の `hence also` —— `Φ^pf` と `Φ^rlf` へ

原文 (FrdI p.115):
> that D is Div-slim [relative to Φ, hence also relative to Φpf, Φrlf].

★★これで `Theorem 6.4, (i)` の Div-slim は**算術で完全に閉じた**。 -/

variable (F Kbar) in
/-- ★★★★★**算術の `𝒟` は `Φ^pf` に関しても Div-slim**（★仮定なし）。 -/
theorem finSubOp_isDivSlim_pf [IsGalois F Kbar] :
    IsDivSlim ((phiGalois F Kbar).pfOn (fun _ => isSharp_effR _)) :=
  isDivSlim_pfOn _ (phiGalois_isDivisorialOn F Kbar) (finSubOp_isDivSlim' F Kbar)

variable (F Kbar) in
/-- ★★★★★**算術の `𝒟` は `Φ^rlf` に関しても Div-slim**（★仮定なし）。

★`hcharInj` は `phiGalois_isCharInj_scMap`（`ArithRlfFlat.lean`、
相対基本等式とトレース写像で作った）である。 -/
theorem finSubOp_isDivSlim_rlf [IsGalois F Kbar] :
    IsDivSlim (phiScOn ℝ≥0 (phiGalois F Kbar)
      (fun {_ _} α => phiGalois_isCharInj_scMap F Kbar α)) :=
  isDivSlim_phiScOn_nnreal _ (fun A => (phiGalois_isPerfFactorial F Kbar A).choose)
    (fun A => (phiGalois_isPerfFactorial F Kbar A).choose_spec)
    (phiGalois_isDivisorialOn F Kbar) (finSubOp_isDivSlim' F Kbar)

end Arith

/-! ### ★出典の紐付け -/

def isDivSlim_of_endo_map_injective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — Φ が自己準同型を区別すれば Div-slim(一般形)",
    sectionId := "frdi-thm-6-4" }

def arithDatum_isDivSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — ArithDatum の水準の Div-slim",
    sectionId := "frdi-thm-6-4" }

def resPlaceOf_eq_of_pullOf_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 引き戻しの一致から素点の制限の一致へ",
    sectionId := "frdi-thm-6-4" }

def finSubOp_isDivSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 𝒟 は Φ に関して Div-slim",
    sectionId := "frdi-thm-6-4" }

def finSubOp_isDivSlim.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "eq_one_of_fixes_prime_of_galoisClosure(数論の側。Chebotarev 不使用)"
      (.inProject "ABC3" "ABC3.Found.NF.eq_one_of_fixes_prime_of_galoisClosure") 115,
    .citation "[ABC3]" "localDeg_eq_one_of_bijective'(自己同型に沿った局所次数は 1)"
      (.inProject "ABC3" "ABC3.Found.Divisor.localDeg_eq_one_of_bijective'") 113,
    .citation "[ABC3]" "isGenSubgroupR_arithDivGroup(各素点に台をもつ正の算術因子)"
      (.inProject "ABC3" "ABC3.Found.Divisor.isGenSubgroupR_arithDivGroup") 113,
    .citation "[ABC3]" "finSub_endo_eq_of_resPlaceOf_eq(仮定 hfaith を落とした形)"
      (.inProject "ABC3" "ABC3.Found.FrdI.finSub_endo_eq_of_resPlaceOf_eq") 115,
    .implicitStep
      ("★訂正: 当初 [Mzk7] Cor 1.1.6 が要ると見積もっていたが、" ++
       "Φ.map e が素点の引き戻しなので圏論の入力は要らない") 115 ]


def finSub_endo_eq_of_resPlaceOf_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 素点の制限が一致する自己射は等しい(仮定なし)",
    sectionId := "frdi-thm-6-4" }

def finSub_endo_eq_of_resPlaceOf_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "eq_one_of_fixes_all_primes(全素イデアルを固定する自己同型は恒等。Chebotarev 不使用)"
      (.inProject "ABC3" "ABC3.Found.NF.eq_one_of_fixes_all_primes") 115,
    .derivation "原文「any automorphism of a number field that fixes all of the valuations … is clearly equal to the identity」" 115 ]

def finSubOp_isDivSlim'.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 𝒟 は Φ に関して Div-slim(仮定なし)",
    sectionId := "frdi-thm-6-4" }

def finSubOp_isDivSlim_pf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 𝒟 は Φ^pf に関して Div-slim(hence also の第 1)",
    sectionId := "frdi-thm-6-4" }

def finSubOp_isDivSlim_rlf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 𝒟 は Φ^rlf に関して Div-slim(hence also の第 2)",
    sectionId := "frdi-thm-6-4" }

def intMapOf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 𝓞 の水準の射(Algebra.id を拾わないよう固定する)",
    sectionId := "frdi-thm-6-4" }

def comap_intMapOf_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 素点の制限の一致は素イデアルの引き戻しの一致",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.FrdI
