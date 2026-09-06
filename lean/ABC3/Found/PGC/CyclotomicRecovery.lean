import ABC3.Skeleton.PGC.Section1
import ABC3.Found.PGC.GaloisTransferContinuous
import Mathlib.RingTheory.RootsOfUnity.AlgebraicallyClosed

/-!
# [pGC] Proposition 1.1 — 円分指標の移送のうち、**体の同型から来る部分**

原文 (pGC p.3):

> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

`Skeleton/PGC/Section1.lean::cyclotomicCharacter_recoverable` は、この主張を
`AssociatedObject.RecoverableFromAbsGal`(§1 冒頭の暗黙の定義)の形

  `∀ K K' (α : Γ_K ≃ₜ* Γ_{K'}), χ_{K'} ∘ α = χ_K`

に写したもので、そこは今も `sorry` である。本ファイルは **`sorry` を持たない**部分だけを
確立する。

## 本ファイルが証明したこと

1. **`cyclotomicCharacter_conj`**(一般形、体でなくても整域でよい)——
   `τ : L ≃+* L'` による共役 `g ↦ τ ∘ g ∘ τ⁻¹` は円分指標を保つ:
   `χ_{L'}(τ ∘ g ∘ τ⁻¹) = χ_L(g)`。mathlib には円分指標の自然性が無いので自前に作る。
2. **`cyclotomicCharacter_galContinuousMulEquiv`**——上を `PAdicLocalField` に落とし、
   体の同型 `α : K ≃ₐ[ℚ_p] K'` が誘導する `Γ_K ≃ₜ* Γ_{K'}`
   (`Found/PGC/GaloisTransferContinuous.lean::galContinuousMulEquiv`)について
   Proposition 1.1 の結論を**無条件に**得る。
3. **`cyclotomicCharacterObject_recoverable_iff`**——`RecoverableFromAbsGal` を
   点ごとの等式 `χ_{K'}(α g) = χ_K(g)` に還元する(`transport` の展開)。
4. **`padicUnits_eq_of_smul_equivariant`**——階数 1 の `ℤ_p`-表現は、**同型類が指標を決める**:
   `ℤ_p(ψ) ≃ ℤ_p(ψ')`(`ℤ_p`-線形かつ同変)ならば `ψ = ψ'`。

## ★4 は原典の「飛躍」の一段を塞ぐ

`Skeleton/PGC/Section1.lean` の docstring は

> その地の文が確立するのは「Γ_K-加群 Z_p(1) の同型類が回復できる」ことであり、
> 本命題が主張するのは「指標 χ が回復できる」こと。この間の一段は原典にない。

と記録し、① 我々のモデル化の誤り / ② 数学が未構築 / ③ 原典側の飛躍 の切り分けを
求めていた。**4 がその一段そのもの**であり、初等的に閉じる(階数 1 なので
`ℤ_p`-線形同型は「`φ 1` 倍」しかなく、整域で約せる)。すなわち
**③ ではない**——原典の論証はこの点では飛んでいない。
残っているのは ② だけ、すなわち「`Γ_K` から `ℤ_p(1)` の同型類を作る」段
(局所 Tate 双対性 = 副有限群の連続 `H²`)である。
`cyclotomicCharacterObject_transport_of_moduleEquiv` がその接続点を明示する。

## 残っている穴(このファイルでは埋めていない)

2 が扱うのは「`α` が体の同型から来る場合」だけである。Proposition 1.1 の内容は
**そうでない `α`**——原典が Introduction で [8] (Jarden–Ritter) を引いて注意する、
`K` と `K'` が同型でないのに `Γ_K ≅ Γ_{K'}` となる場合——にある。そこは
`cyclotomicCharacterObject_transport_of_moduleEquiv` の仮説
(`ℤ_p(χ_K)` と `ℤ_p(χ_{K'} ∘ α)` が同変同型)が未証明のまま残る。

## ★2026-09-06 の在庫測定(次のノードの候補)

* **`μ_{p^n}` と Lubin–Tate の `ψ`-捩れ点を同定する橋は、この木に無い。**
  `Found/PGC/LubinTateDistinguishedSeparable.lean::iteratedLubinTatePsiTorsionPoints`
  は一般の `K`・一般の素元 `π`・一般の Lubin–Tate 級数について立てられており、
  `K = ℚ_p`・`f = (1+X)^p − 1` に特化して `x ↦ x + 1` で `μ_{p^n}` と同定する補題は
  `.cache/decl-index.txt` に 1 件も無い(`iteratedLubinTatePsiTorsionPoints` を含む宣言 14 件を
  すべて見たが、`rootsOfUnity`/`IsPrimitiveRoot` に触れるものは無い)。
  `Found/PGC/QpRootsOfUnity.lean` は `ℓ ≠ p` の**素な**乗根の話で、`p` 冪乗根は扱っていない。
* **仮に橋を架けて `Gal(K(μ_{p^n})/K) ≃ (ℤ/p^n)ˣ` を得ても、それだけでは
  Proposition 1.1 は出ない。** `RecoverableFromAbsGal` は
  `cyclotomicCharacterObject_recoverable_iff` が示すとおり
  「**任意の**位相群同型 `α` について `χ_{K'}(α g) = χ_K(g)`」であり、
  各 `K` ごとの `Gal(K(μ_{p^n})/K)` の**構造**を知っても、抽象群の同型 `α` の振る舞いは
  制約されない。必要なのは `Γ_K` の中に `χ` を**群論的に特徴づける**手続き
  (円分子 `Λ(Γ_K) := Hom(H²(Γ_K, ℤ/p^n), ℤ/p^n) ≅ μ_{p^n}`)である。
* **`ker(χ mod p^n)` だけなら経路 C 流の計数で届く見込みがある**が、それでは足りない。
  `Found/PGC/KummerDuality.lean::card_contHom_eq_index_powRange` は `μ_n ⊆ F` の体で
  `N_n(Γ_F) = [Fˣ : (Fˣ)ⁿ]` を与えるので、`n = p^m` でも
  「`μ_{p^m} ⊆ L_H` かどうか」を `N_{p^m}(H)` と `[L_H:ℚ_p]` から読める可能性がある
  (`Found/PGC/PrincipalUnitsRank.lean`・`UnitsSplit.lean` が `Lˣ` の構造を持っている)。
  これで `ker(χ_m)` は回復できる。★しかしそこで止まる——核だけを知っても
  `χ_m` は像 `⊆ (ℤ/p^m)ˣ` の自己同型のぶんだけ不定であり、`g ↦ ord(χ_m(g))` を
  すべての `m`・`g` について知っても同じ不定性が残る。
  この不定性を切るのが円分子の**標準的な生成元**、すなわち `H²` の入力である。

## 逸脱

無し。1〜4 はいずれも原典の主張の**部分**であるか、原典が暗黙にしていた一段の明示である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Meta

/-! ## 1. 円分指標は共役で保たれる(一般形) -/

/-- **円分指標の自然性**——`τ : L ≃+* L'` による共役は `χ` を変えない。

`χ_{L'}(τ ∘ g ∘ τ⁻¹) = χ_L(g)`。

mathlib の `cyclotomicCharacter` にはこの自然性が無い(2026-09-06 実測)。
証明は `p^n` 乗根に対する仕様 `cyclotomicCharacter.spec` を両側で使うだけ:
`L'` の原始 `p^n` 乗根 `ζ` を取ると、左辺は `ζ ↦ ζ^a`、右辺は
`τ⁻¹ζ ↦ (τ⁻¹ζ)^b` を `τ` で押し出して `ζ ↦ ζ^b` を与えるので、
`ζ^a = ζ^b` から `a = b`(`IsPrimitiveRoot.pow_inj`)。 -/
theorem cyclotomicCharacter_conj {L L' : Type*} [CommRing L] [IsDomain L]
    [CommRing L'] [IsDomain L'] (p : ℕ) [hp : Fact p.Prime]
    [∀ i, HasEnoughRootsOfUnity L (p ^ i)] [∀ i, HasEnoughRootsOfUnity L' (p ^ i)]
    (τ : L ≃+* L') (g : L ≃+* L) :
    cyclotomicCharacter L' p (τ.symm.trans (g.trans τ)) = cyclotomicCharacter L p g := by
  apply Units.ext
  apply PadicInt.ext_of_toZModPow.mp
  intro n
  haveI : NeZero (p ^ n) := ⟨pow_ne_zero _ hp.out.ne_zero⟩
  obtain ⟨ζ, hζ⟩ := HasEnoughRootsOfUnity.exists_primitiveRoot L' (p ^ n)
  set a := (PadicInt.toZModPow n) ((cyclotomicCharacter L' p (τ.symm.trans (g.trans τ))) : ℤ_[p])
    with ha
  set b := (PadicInt.toZModPow n) ((cyclotomicCharacter L p g) : ℤ_[p]) with hb
  apply ZMod.val_injective
  refine hζ.pow_inj (ZMod.val_lt a) (ZMod.val_lt b) ?_
  have h1 : (τ.symm.trans (g.trans τ)) ζ = ζ ^ a.val :=
    cyclotomicCharacter.spec (n := n) p _ ζ hζ.pow_eq_one
  have hs : (τ.symm ζ) ^ p ^ n = 1 := by
    rw [← map_pow, hζ.pow_eq_one, map_one]
  have h2 : g (τ.symm ζ) = (τ.symm ζ) ^ b.val :=
    cyclotomicCharacter.spec (n := n) p g (τ.symm ζ) hs
  have h3 : (τ.symm.trans (g.trans τ)) ζ = τ (g (τ.symm ζ)) := rfl
  rw [h1] at h3
  rw [h3, h2, map_pow, RingEquiv.apply_symm_apply]

/-! ## 2. 体の同型から来る `α` については Proposition 1.1 が成り立つ -/

variable {p : ℕ} [Fact p.Prime]

/-- **体の同型が誘導する `Γ_K ≃ₜ* Γ_{K'}` は円分指標を保つ**(延長 `ᾱ` を指定する版)。

`galContinuousMulEquivOf α ᾱ hfwd g` は定義上 `ᾱ ∘ g ∘ ᾱ⁻¹` なので、
`cyclotomicCharacter_conj` をそのまま当てればよい。 -/
theorem cyclotomicCharacter_galContinuousMulEquivOf {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) (g : K.absGal) :
    cyclotomicCharacter K'.closure p (galContinuousMulEquivOf α ᾱ hfwd g).toRingEquiv
      = cyclotomicCharacter K.closure p g.toRingEquiv :=
  cyclotomicCharacter_conj p ᾱ g.toRingEquiv

/-- **体の同型が誘導する `Γ_K ≃ₜ* Γ_{K'}` は円分指標を保つ**(既定の延長)。 -/
theorem cyclotomicCharacter_galContinuousMulEquiv {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (g : K.absGal) :
    cyclotomicCharacter K'.closure p (galContinuousMulEquiv α g).toRingEquiv
      = cyclotomicCharacter K.closure p g.toRingEquiv :=
  cyclotomicCharacter_conj p (extendToClosure α) g.toRingEquiv

/-! ## 3. `RecoverableFromAbsGal` の展開 -/

/-- `cyclotomicCharacterObject` の `transport` を展開した形。

`transport α f = f ∘ α⁻¹` なので、`transport α (obj K) = obj K'` は
「`α` に沿って χ が対応する」という点ごとの等式と同値。 -/
theorem cyclotomicCharacterObject_transport_eq_iff {K K' : PAdicLocalField p}
    (α : ContinuousMulEquiv K.absGal K'.absGal) :
    (cyclotomicCharacterObject (p := p)).transport α
        ((cyclotomicCharacterObject (p := p)).obj K)
      = (cyclotomicCharacterObject (p := p)).obj K'
    ↔ ∀ g : K.absGal, cyclotomicCharacter K'.closure p (α g).toRingEquiv
        = cyclotomicCharacter K.closure p g.toRingEquiv := by
  constructor
  · intro h g
    have h2 : (cyclotomicCharacterObject (p := p)).transport α
        ((cyclotomicCharacterObject (p := p)).obj K) (α g)
        = (cyclotomicCharacterObject (p := p)).obj K' (α g) := congrFun h (α g)
    have h3 : (cyclotomicCharacterObject (p := p)).transport α
        ((cyclotomicCharacterObject (p := p)).obj K) (α g)
        = cyclotomicCharacter K.closure p g.toRingEquiv := by
      show cyclotomicCharacter K.closure p (α.toMulEquiv.symm (α g)).toRingEquiv = _
      simp
    exact (h3.symm.trans h2).symm
  · intro h
    funext g'
    show cyclotomicCharacter K.closure p (α.toMulEquiv.symm g').toRingEquiv
      = cyclotomicCharacter K'.closure p g'.toRingEquiv
    have := h (α.toMulEquiv.symm g')
    rw [show α (α.toMulEquiv.symm g') = g' by simp] at this
    exact this.symm

/-- **Proposition 1.1 の点ごとの形**——`RecoverableFromAbsGal` は
「どんな `α` でも `χ_{K'}(α g) = χ_K(g)`」と同値。 -/
theorem cyclotomicCharacterObject_recoverable_iff :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal
    ↔ ∀ (K K' : PAdicLocalField p) (α : ContinuousMulEquiv K.absGal K'.absGal)
        (g : K.absGal),
        cyclotomicCharacter K'.closure p (α g).toRingEquiv
          = cyclotomicCharacter K.closure p g.toRingEquiv := by
  constructor
  · intro h K K' α
    exact (cyclotomicCharacterObject_transport_eq_iff α).mp (h α)
  · intro h K K' α
    exact (cyclotomicCharacterObject_transport_eq_iff α).mpr (h K K' α)

/-- **体の同型から来る `α` については Proposition 1.1 が成り立つ**(対象の言葉で)。

Proposition 1.2 の退化(`Check/PGC/Prop12Degenerate.lean`)を作ったのは
まさにこの経路の `α` だった。Proposition 1.1 はその経路では退化しない。 -/
theorem cyclotomicCharacterObject_transport_galContinuousMulEquiv
    {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    (cyclotomicCharacterObject (p := p)).transport (galContinuousMulEquiv α)
        ((cyclotomicCharacterObject (p := p)).obj K)
      = (cyclotomicCharacterObject (p := p)).obj K' :=
  (cyclotomicCharacterObject_transport_eq_iff (galContinuousMulEquiv α)).mpr
    (fun g => cyclotomicCharacter_galContinuousMulEquiv α g)

/-! ## 4. 階数 1 では「同型類」が「指標」を決める -/

/-- **階数 1 の `ℤ_p`-表現は同型類で指標が決まる。**

`ψ, ψ' : G → ℤ_[p]ˣ` について、`ℤ_p` 上の線形同型 `φ` が
`φ (ψ g • x) = ψ' g • φ x` を満たせば `ψ = ψ'`。

証明: `x = 1` に入れると `ψ g * φ 1 = ψ' g * φ 1`。`φ` は同型なので `φ 1 ≠ 0`、
`ℤ_p` は整域なので約せる。

★これが `Skeleton/PGC/Section1.lean` の docstring が
「原典にない一段」と記録していた
「`ℤ_p(1)` の同型類が回復できる ⟹ 指標 χ が回復できる」の一段である。
階数 1 なので初等的に閉じる。 -/
theorem padicUnits_eq_of_smul_equivariant {G : Type*} {ψ ψ' : G → ℤ_[p]ˣ}
    (φ : ℤ_[p] ≃ₗ[ℤ_[p]] ℤ_[p])
    (h : ∀ (g : G) (x : ℤ_[p]), φ ((ψ g : ℤ_[p]) • x) = (ψ' g : ℤ_[p]) • φ x) :
    ψ = ψ' := by
  have hφ1 : φ 1 ≠ 0 := fun h0 => one_ne_zero (φ.map_eq_zero_iff.mp h0)
  funext g
  apply Units.ext
  have h1 : ((ψ g : ℤ_[p]) * φ 1) = ((ψ' g : ℤ_[p]) * φ 1) := by
    have := h g 1
    rwa [map_smul, smul_eq_mul, smul_eq_mul] at this
  exact mul_right_cancel₀ hφ1 h1

/-- **接続点**——「`ℤ_p(1)` の同型類が `α` に沿って対応する」から
Proposition 1.1 の結論が出る。

仮説は「`Γ_K` の作用を `α` で移した `ℤ_p(χ_{K'} ∘ α)` が `ℤ_p(χ_K)` と
`ℤ_p`-線形同変同型である」こと——これが局所 Tate 双対性が与えるはずのもので、
**まだ証明されていない**(mathlib にも本リポジトリにも副有限群の連続 `H²` の
道具が無い。`Skeleton/PGC/Section1.lean::cyclotomicCharacter_recoverable.needs` 参照)。
本補題はその一段だけを仮定に切り出し、残り(同型類 ⟹ 指標)を消す。 -/
theorem cyclotomicCharacterObject_transport_of_moduleEquiv
    (h : ∀ (K K' : PAdicLocalField p) (α : ContinuousMulEquiv K.absGal K'.absGal),
      ∃ φ : ℤ_[p] ≃ₗ[ℤ_[p]] ℤ_[p], ∀ (g : K.absGal) (x : ℤ_[p]),
        φ ((cyclotomicCharacter K.closure p g.toRingEquiv : ℤ_[p]) • x)
          = (cyclotomicCharacter K'.closure p (α g).toRingEquiv : ℤ_[p]) • φ x) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal := by
  refine cyclotomicCharacterObject_recoverable_iff.mpr ?_
  intro K K' α g
  obtain ⟨φ, hφ⟩ := h K K' α
  have := padicUnits_eq_of_smul_equivariant (G := K.absGal)
    (ψ := fun g => cyclotomicCharacter K.closure p g.toRingEquiv)
    (ψ' := fun g => cyclotomicCharacter K'.closure p (α g).toRingEquiv) φ hφ
  exact (congrFun this g).symm

end ABC3.Found.PGC
