/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Remark151Chain

/-!
# [GenEll] Remark 1.5.1 —— **`ratTower n` は `ℤ` の局所化である**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★これで点の持ち上げが作れる

`Remark151Chain.lean` は「点が段 `n` のモデルへ持ち上がる」ことを**仮定**として受けていた。
★その仮定を消すには `ratTower n ⟶ 𝓞_F[1/n!]` が要る。

`ratTower n = Subring.closure {(n!)⁻¹} ⊆ ℚ` は**定義が閉包**なので、
そのままでは普遍性が使えない。★★そこで

    IsLocalization (Submonoid.powers ((n! : ℕ) : ℤ)) (ratTower n)

を示す。★★★これが付けば `IsLocalization.lift` で
「`n!` が可逆な任意の `ℤ`-代数」への環準同型が一斉に出る。

## ★★★★★★機構

* `map_units` —— `(n!)⁻¹` が閉包に入っているので `n!` は可逆。冪も同様。
* `surj` —— ★**`ratTower n` の元は `a / (n!)^k` の形をしている**。
  これを `Subring.closure_induction` で示す（`ratTower_exists_num`）。
* `exists_of_eq` —— `ℤ ↪ ℚ` なので単射。`c = 1` で済む。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★★★★★`ratTower n` の元の形 -/

/-- ★★★★★★**`ratTower n` の元は `a / (n!)^k` の形をしている**。

★`Subring.closure_induction` で 6 通り（生成元・0・1・和・逆元・積）を潰す。
★★指数は和のところで `k + l` に取るのが要点である。 -/
theorem ratTower_exists_num (n : ℕ) (z : ratTower n) :
    ∃ (a : ℤ) (k : ℕ), (z : ℚ) * ((Nat.factorial n : ℕ) : ℚ) ^ k = (a : ℚ) := by
  have hz : (z : ℚ) ∈ ratTower n := z.2
  refine Subring.closure_induction (p := fun y _ =>
    ∃ (a : ℤ) (k : ℕ), y * ((Nat.factorial n : ℕ) : ℚ) ^ k = (a : ℚ)) ?_ ?_ ?_ ?_ ?_ ?_ hz
  · rintro x rfl
    refine ⟨1, 1, ?_⟩
    have hne : ((Nat.factorial n : ℕ) : ℚ) ≠ 0 := by
      exact_mod_cast Nat.cast_ne_zero.2 (Nat.factorial_ne_zero n)
    rw [pow_one, inv_mul_cancel₀ hne]
    norm_num
  · exact ⟨0, 0, by simp⟩
  · exact ⟨1, 0, by simp⟩
  · rintro x y _ _ ⟨a, k, ha⟩ ⟨b, l, hb⟩
    refine ⟨a * ((Nat.factorial n : ℕ) : ℤ) ^ l + b * ((Nat.factorial n : ℕ) : ℤ) ^ k, k + l, ?_⟩
    push_cast
    rw [← ha, ← hb]
    ring
  · rintro x _ ⟨a, k, ha⟩
    exact ⟨-a, k, by push_cast; rw [← ha]; ring⟩
  · rintro x y _ _ ⟨a, k, ha⟩ ⟨b, l, hb⟩
    refine ⟨a * b, k + l, ?_⟩
    push_cast
    rw [← ha, ← hb]
    ring

/-! ## ★★★★`n!` は `ratTower n` で可逆 -/

/-- ★`(n!)⁻¹` を `ratTower n` の元として。 -/
noncomputable def invFac (n : ℕ) : ratTower n :=
  ⟨((Nat.factorial n : ℕ) : ℚ)⁻¹, Subring.subset_closure rfl⟩

theorem fac_mul_invFac (n : ℕ) :
    algebraMap ℤ (ratTower n) ((Nat.factorial n : ℕ) : ℤ) * invFac n = 1 := by
  have hne : ((Nat.factorial n : ℕ) : ℚ) ≠ 0 := by
    exact_mod_cast Nat.cast_ne_zero.2 (Nat.factorial_ne_zero n)
  apply Subtype.ext
  push_cast [invFac]
  rw [mul_inv_cancel₀ hne]

theorem isUnit_fac (n : ℕ) :
    IsUnit (algebraMap ℤ (ratTower n) ((Nat.factorial n : ℕ) : ℤ)) :=
  IsUnit.of_mul_eq_one _ (fac_mul_invFac n)

theorem isUnit_fac_pow (n k : ℕ) :
    IsUnit (algebraMap ℤ (ratTower n) (((Nat.factorial n : ℕ) : ℤ) ^ k)) := by
  rw [map_pow]
  exact (isUnit_fac n).pow k

/-! ## ★★★★★★★★★★局所化であること -/

/-- ★★★★★★★★★★**`ratTower n` は `ℤ` の `n!` による局所化である**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★これで `IsLocalization.lift` が使えるようになり、
「`n!` が可逆な任意の `ℤ`-代数」への環準同型が一斉に出る。 -/
instance ratTower_isLocalization (n : ℕ) :
    IsLocalization (Submonoid.powers ((Nat.factorial n : ℕ) : ℤ)) (ratTower n) where
  map_units := by
    rintro ⟨y, k, rfl⟩
    exact isUnit_fac_pow n k
  surj := by
    intro z
    obtain ⟨a, k, ha⟩ := ratTower_exists_num n z
    refine ⟨⟨a, ⟨((Nat.factorial n : ℕ) : ℤ) ^ k, k, rfl⟩⟩, ?_⟩
    apply Subtype.ext
    push_cast
    exact ha
  exists_of_eq := by
    intro x y h
    refine ⟨1, ?_⟩
    have hq : ((x : ℚ)) = ((y : ℚ)) := by
      have := congrArg (fun t : ratTower n => (t : ℚ)) h
      simpa using this
    have hxy : x = y := by exact_mod_cast hq
    rw [hxy]

/-- ★★★★★★**`n!` が可逆な `ℤ`-代数への環準同型**（普遍性）。 -/
noncomputable def ratTowerLift (n : ℕ) (B : Type) [CommRing B]
    (hinv : IsUnit (algebraMap ℤ B ((Nat.factorial n : ℕ) : ℤ))) : ratTower n →+* B :=
  IsLocalization.lift (M := Submonoid.powers ((Nat.factorial n : ℕ) : ℤ))
    (g := algebraMap ℤ B) (by rintro ⟨y, k, rfl⟩; rw [map_pow]; exact hinv.pow k)

/-! ## ★★★★★★★★★★点の持ち上げ -/

variable (F : Type) [Field F] [NumberField F]

/-- ★★★★★★★★**点を段 `n` のモデルへ持ち上げる**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`Spec 𝓞_F[1/n!] ⟶ X ×_ℤ ℤ[1/n!]` を作る。

★成分は 2 つ: `Spec 𝓞_F[1/n!] ⟶ Spec ℤ[1/n!]`（`ratTowerLift` の `Spec`）と
`Spec 𝓞_F[1/n!] ⟶ Spec 𝓞_F ⟶ X`。
★★両立は `Spec ℤ` が終対象なので自動。 -/
noncomputable def liftPointToBc {n : ℕᵒᵖ} (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    (xF : specRingOfIntegers F ⟶ X) :
    Spec (CommRingCat.of A) ⟶ bcObj f n :=
  pullback.lift
    (Spec.map (CommRingCat.ofHom (ratTowerLift n.unop A hinv)))
    (Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF)
    (specZIsTerminal.hom_ext _ _)

@[simp] theorem liftPointToBc_snd {n : ℕᵒᵖ} (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    (xF : specRingOfIntegers F ⟶ X) :
    liftPointToBc F A hinv f xF ≫ pullback.snd (overRatTowerDiagram.obj n).hom f
      = Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF :=
  pullback.lift_snd _ _ _

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— 持ち上げの仮定が消えた形**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`Remark151Chain.lean` の `conductorADiv_fin_eq_of_lift` は
「点が段 `n` へ持ち上がる」を仮定として受けていた。
★★**`ratTower n` が局所化であることが取れたので、その持ち上げは構成できる**
——本定理では仮定が `hinv`（`n!` が `A` で可逆）だけになる。

★★★残るのは `xF'`（`X'` 側の点）がその持ち上げと両立するという仮定 `hlift'` である
——これは原文では `X_ℚ ≅ X′_ℚ` が点の対応を与える段に当たり、
`𝓞_F` へ戻すのに `X'` の固有性（付値判定法）を使う。 -/
theorem conductorADiv_fin_eq_of_isUnit
    (M : Submonoid (NumberField.RingOfIntegers F)) (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A] [IsLocalization M A]
    (v : FinitePlace F) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    {n : ℕᵒᵖ}
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ)
    (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X')
    (hlift' : Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF'
      = liftPointToBc F A hinv f xF ≫ φ ≫
        pullback.snd (overRatTowerDiagram.obj n).hom f')
    (hI : pullbackIdeal F iZ.ker xF ≠ 0) (hJ : pullbackIdeal F iZ'.ker xF' ≠ 0) :
    (conductorADiv F iZ.ker xF).fin v = (conductorADiv F iZ'.ker xF').fin v :=
  conductorADiv_fin_eq_of_lift F M A v hM f f' iZ iZ' φ ψ hsq
    (liftPointToBc F A hinv f xF) xF xF'
    (liftPointToBc_snd F A hinv f xF).symm hlift' hI hJ

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` はまだ置かない。** 残っているのは
`hlift'`——`X'` 側の点が持ち上げと両立すること——であり、
原文では `X_ℚ ≅ X′_ℚ` が点の対応を与える段に当たる。
`𝓞_F` へ戻すのに `X'` の固有性（付値判定法）を使う。 -/

def ratTower_isLocalization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(ratTower n = ℤ[1/n!] が ℤ の局所化であること)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_fin_eq_of_isUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(hagree——点の持ち上げを構成した形。X′ 側の点の対応は仮定)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_fin_eq_of_isUnit.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "conductorADiv_fin_eq_of_lift(鎖の全体)"
      (.inProject "ABC3" "ABC3.Found.GenEll.conductorADiv_fin_eq_of_lift") 9,
    .citation "[mathlib]" "IsLocalization.lift(局所化の普遍性)"
      (.inMathlib "IsLocalization.lift") 9,
    .citation "[mathlib]" "Subring.closure_induction(閉包の元の形を決める)"
      (.inMathlib "Subring.closure_induction") 9,
    .implicitStep
      ("★ratTower n = Subring.closure {(n!)⁻¹} は定義が閉包なので、" ++
       "そのままでは普遍性が使えない。★★元が a/(n!)^k の形をしていることを" ++
       "closure_induction で示して IsLocalization を立てた") 9,
    .implicitStep
      ("★★★残る仮定 hlift': X′ 側の点が持ち上げと両立すること。" ++
       "★原文では X_ℚ ≅ X′_ℚ が点の対応を与える段に当たり、" ++
       "𝓞_F へ戻すのに X′ の固有性(付値判定法)を使う。" ++
       "★★mathlib には固有射の付値判定法があるので、欠落ではなく我々の作業である") 9 ]

end ABC3.Found.GenEll
