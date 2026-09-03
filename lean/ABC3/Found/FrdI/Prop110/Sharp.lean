/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.GroupTheory.MonoidLocalization.GrothendieckGroup
import Mathlib.NumberTheory.PrimeCounting
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.MonoidPrime
import ABC3.Found.FrdI.Prop110.PreStep

/-!
# Prop110 —— sharp モノイド・(iv) の逆向き・(vi) の前半

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ### ★sharp モノイドの 2 つの補題 —— (iv) の逆向きに要る

★**どちらも「和が 0 なら可逆、可逆なら 0」という sharp の一行の帰結**だが、
**使う形が違う**ので別々に置く。

★**原文はこれらを一度も書かない。** `Definition 1.1, (i)` で `divisorial` に
`sharp` を入れておいて、以後は暗黙に使う。★**「sharp を仮定に入れる」という
一回の判断が、下流の何箇所で効いているかは、形式化するまで見えない。**
-/

/-- ★**sharp なら、和が 0 の両項は 0**。 -/
theorem eq_zero_of_add_eq_zero_of_isSharp {M : Type w} [AddCommMonoid M]
    (hs : IsSharp M) {a b : M} (h : a + b = 0) : a = 0 ∧ b = 0 :=
  ⟨hs a ⟨⟨a, b, h, by rw [add_comm]; exact h⟩, rfl⟩,
   hs b ⟨⟨b, a, by rw [add_comm]; exact h, h⟩, rfl⟩⟩

/-- ★**sharp なら、`n • a = 0`(`n ≥ 1`)から `a = 0`**。

★`n • a = a + (n-1) • a` なので `a` は可逆。 -/
theorem eq_zero_of_nsmul_eq_zero_of_isSharp {M : Type w} [AddCommMonoid M]
    (hs : IsSharp M) {n : ℕ} (hn : 0 < n) {a : M} (h : n • a = 0) : a = 0 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [succ_nsmul] at h
  exact (eq_zero_of_add_eq_zero_of_isSharp hs h).2

/-! ### ★(iv) の逆向き —— prime-Frobenius ⟹ irreducible

★**任意の分解 `φ = β ≫ α` に対して片方が同型**を示す。分解は**任意の射**なので、
`β`・`α` が Frobenius 型であることは**与えられていない**。そこで 3 つを別々に出す:

1. **base-isomorphism**: `Base φ = Base β ≫ Base α` が同型で `𝒟` は totally epimorphic
   なので、**両方が同型**(`isIso_of_isIso_comp`)
2. **isometric**: `Div φ = 0`(Frobenius 型)を `Div_comp` で開くと
   `Φ.map (Base β) (Div α) + (degFr α) • Div β = 0`。
   ★**sharp なので両項が 0**、さらに `Φ.map` は**常に単射**なので `Div α = 0`、
   `n • x = 0` から `Div β = 0`。
3. **co-angular**: `𝒞` が isotropic 型なので `Proposition 1.4, (i)` から自動

★あとは次数。`degFr α * degFr β = degFr φ` が素数なので
`pnat_irreducible_iff_prime` の `⟸` 向きから **片方が 1**。
次数 1 ＋ LB-invertible ＋ base-iso なら `Proposition 1.4, (iii)` で同型。

★★**ここで sharp が 2 回効いている**(和の分解と `n` 倍の消去)。
★**原文は (iv) を「Proposition 1.7, (v) と ℕ≥1 の well-known structure から
immediately」と書くだけで、この 3 段は書いていない。**
-/

include P in
/-- **(iv) の逆向き** —— prime-Frobenius なら irreducible。

★`𝒞` が isotropic 型であることを使う(原文の `Proposition 1.4, (i)` 経由)。 -/
theorem prop_1_10_iv_mp (F : FrobenioidCore P) {A B : C} (hA : IsIsotropic P A)
    (φ : A ⟶ B) (hp : IsPrimeFrobenius P φ) : IsIrreducibleMor φ := by
  refine ⟨fun h => ?_, fun X β α hfac => ?_⟩
  · haveI := h
    have h2 := hp.2
    rw [degFr_of_isIso P φ] at h2
    exact Nat.not_prime_one (by simpa using h2)
  · -- base-isomorphism を両方に降ろす
    have hb : IsIso (P.Base β ≫ P.Base α) := by
      rw [← P.Base_comp, ← hfac]; exact hp.1.2
    haveI := hb
    obtain ⟨hbβ, hbα⟩ := isIso_of_isIso_comp P.totEpiD (P.Base β) (P.Base α)
    -- isometric を両方に降ろす
    have hsum : Φ.map (P.Base β) (P.Div α) + (P.degFr α : ℕ) • P.Div β = 0 := by
      rw [← P.Div_comp, ← hfac]; exact hp.1.1.2
    obtain ⟨e1, e2⟩ := eq_zero_of_add_eq_zero_of_isSharp (P.divisorial _).2 hsum
    have hdβ : P.Div β = 0 :=
      eq_zero_of_nsmul_eq_zero_of_isSharp (P.divisorial _).2 (P.degFr α).pos e2
    have hdα : P.Div α = 0 :=
      Φ.map_injective (P.Base β) (by rw [e1, map_zero])
    -- 次数の分解
    have hd : P.degFr φ = P.degFr α * P.degFr β := by rw [hfac, P.degFr_comp]
    rcases ((pnat_irreducible_iff_prime (P.degFr φ)).mpr hp.2).2 _ _ hd with h | h
    · -- degFr α = 1 ⟹ α は同型
      -- ★`α : X ⟶ B` の域 `X` は `β : A ⟶ X` を通じて `A` から到達するので isotropic
      exact Or.inr (prop_1_4_iii P F α
        ⟨prop_1_4_i P α (fun _ f => F.isotropicClosed f (F.isotropicClosed β hA)), hdα⟩
        ⟨h, hbα⟩)
    · -- degFr β = 1 ⟹ β は同型
      exact Or.inl (prop_1_4_iii P F β
        ⟨prop_1_4_i P β (fun _ f => F.isotropicClosed f hA), hdβ⟩ ⟨h, hbβ⟩)

include P in
/-- ★★**`Proposition 1.10, (iv)` の前半の完成形**(iff)。

★★**仮定は「域 `A` が isotropic」だけ**である。

原文 (FrdI p.34):
> (iv) A morphism of Frobenius type with isotropic domain is a prime-Frobe-

★**監査以前は `∀ X : C, IsIsotropic P X`(圏全体)を仮定していた** ——
原文が要求しない仮定を足していた。`F.isotropicClosed` で域から伝播させれば足りる。 -/
theorem prop_1_10_iv (F : FrobenioidCore P) {A B : C} (hA : IsIsotropic P A)
    (φ : A ⟶ B) (hφ : IsFrobeniusType P φ) :
    IsPrimeFrobenius P φ ↔ IsIrreducibleMor φ :=
  ⟨prop_1_10_iv_mp P F hA φ,
   prop_1_10_iv_mpr P F (fun _ f => F.isotropicClosed f hA) φ hφ⟩

/-! ### ★(vi) の前半の核心 —— base-identity 自己射を pre-step に沿って降ろす

原文 (FrdI p.36、目視):
> since φC is a base-identity endomorphism of Frobenius degree d, and γ is a pre-
> step, it follows [cf. Remark 1.1.1] that φB is also a base-identity endomorphism of
> Frobenius degree d.

★**この一段だけを取り出す**。四角形 `γ ≫ φ_B = φ_C ≫ γ` があり、
`φ_C` が base-identity で次数 `d`、`γ` が pre-step のとき、
★**`φ_B` も base-identity で次数 `d`**。

★**証明の中身**:
- **次数**: `degFr (γ ≫ φ_B) = degFr φ_B * degFr γ = degFr φ_B`(`γ` は pre-step ⟹ 次数 1)、
  `degFr (φ_C ≫ γ) = degFr γ * degFr φ_C = d`。よって `degFr φ_B = d`。
- **base-identity**: `Base γ ≫ Base φ_B = Base φ_C ≫ Base γ = Base γ`(`φ_C` が base-identity)。
  ★**`𝒟` が totally epimorphic なので `Base γ` は epi**、よって `Base φ_B = 𝟙`。

★★**原文の「[cf. Remark 1.1.1]」は次数の部分だけを指している。**
base-identity の部分は **`𝒟` の totally epimorphic 性**から出るもので、
`Remark 1.1.1`(合成則)ではない。★**1 つの括弧が 2 つの理由を覆っていた。**
-/

include P in
/-- **(vi) の前半の核心** —— base-identity 自己射を pre-step に沿って降ろす。

★次数は `Remark 1.1.1`(合成則)から、base-identity は **`𝒟` の totally epimorphic 性**から。
★**原文の括弧 `[cf. Remark 1.1.1]` は前者だけを指している。** -/
theorem prop_1_10_vi_descend {B Cc : C} (γ : B ⟶ Cc) (hγ : IsPreStep P γ)
    (φC : Cc ⟶ Cc) (φB : B ⟶ B) (hbC : IsBaseIdentity P φC)
    {d : ℕ+} (hdC : P.degFr φC = d)
    (hsq : φB ≫ γ = γ ≫ φC) :
    IsBaseIdentity P φB ∧ P.degFr φB = d := by
  haveI : IsIso (P.Base γ) := hγ.2
  constructor
  · -- base-identity: ★`γ` が pre-step ⟹ `Base γ` は**同型**なので右から消せる
    have h := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp, show P.Base φC = P.Base (𝟙 Cc) from hbC,
      P.Base_id, Category.comp_id] at h
    show P.Base φB = P.Base (𝟙 B)
    rw [P.Base_id]
    exact (cancel_mono (P.Base γ)).mp (by rw [h, Category.id_comp])
  · -- 次数: 合成則から
    have h := congrArg P.degFr hsq
    rw [P.degFr_comp, P.degFr_comp, hγ.1, hdC, one_mul, mul_one] at h
    exact h

/-! ## ★★`Definition 1.3, (iii), (d)` の圏同値を**消費する**

★**測定(2026-08-16)**: `coaPreUnderEquiv` / `coaPreOverEquiv` は
`Prop15`(`𝔽_Φ` の場合)と `Prop19`(`𝒞^istr` への移送)で**構成**されているが、
★**下流で一度も消費されていない**(`grep` で確認)。

★**`Proposition 1.10` の残り 2 件((iii) の後半、(vi) の前半)は
どちらもこの圏同値を待っている。** だから**橋を 1 本架ける**。

★**橋の形**: `(𝒞^coa-pre)_{Cc} ⥤ Order(Φ(Cc))^opp` が**充満**であることから、
`Div` の順序関係 `γ の不変量 ≤ γ′ の不変量` を、
★**射の因子分解 `β ≫ γ = γ′`** に変換する。

★**なぜ `op` が要るか**: 関手は反変(`Order(...)^opp` へ行く)。
`F(Z_{γ′}) ⟶ F(Z_γ)` は `op γ′_inv ⟶ op γ_inv`、すなわち
`γ_inv ⟶ γ′_inv` の `op`、すなわち **`γ_inv ≤ γ′_inv`**。
★**「大きいほうから小さいほうへ射がある」の向きが、`op` で 1 回ひっくり返る。**
-/

include P in
/-- ★★**`Definition 1.3, (iii), (d)` の第2の圏同値の初の消費**。

co-angular pre-step の `Div` 不変量に順序関係があれば、**因子分解が存在する**。

★原文(p.36)の
> the portion of assertion (ii) concerning the relationship between Div(γ), Div(γ′)
> implies, in light of the second equivalence of categories of Definition 1.3, (iii), (d),
> that γ′ factors through γ
の「in light of the second equivalence」の実体である。

★**`[IsIso (P.Base γ)]` を instance binder で取っている理由**: `inv (P.Base γ)` が
**文の中**に現れるので、`hγs.2` を証明の中で `haveI` しても遅い。
★**「仮定に含まれている事実を、インスタンスとしても渡す」必要がある** ——
`IsPreStep` は構造体の連言であってインスタンスではないから。
★これは分類表 #5(インスタンス合成)の**設計側の現れ**である。 -/
theorem coaPre_factor_of_mle (G : Frobenioid P) {Cc B B' : C}
    (γ : B ⟶ Cc) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (γ' : B' ⟶ Cc) (hγ'c : IsCoAngular P γ') (hγ's : IsPreStep P γ')
    [IsIso (P.Base γ)] [IsIso (P.Base γ')]
    (hle : MLe (Φ.map (inv (P.Base γ)) (P.Div γ) : Φ.val (P.toElem.obj Cc).base)
      (Φ.map (inv (P.Base γ')) (P.Div γ'))) :
    ∃ β : B' ⟶ B, IsCoAngular P β ∧ IsPreStep P β ∧ β ≫ γ = γ' := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv Cc
  set Zγ : Over (⟨Cc⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨Cc⟩ from ⟨γ, hγc, hγs⟩) with hZγ
  set Zγ' : Over (⟨Cc⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B'⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨Cc⟩ from ⟨γ', hγ'c, hγ's⟩) with hZγ'
  have hmor : (coaPreOverFunctor P Cc).obj Zγ' ⟶ (coaPreOverFunctor P Cc).obj Zγ :=
    (homOfLE (show (toOrderCat (Φ.map (inv (P.Base γ)) (P.Div γ)))
      ≤ toOrderCat (Φ.map (inv (P.Base γ')) (P.Div γ')) from hle)).op
  let f := (coaPreOverFunctor P Cc).preimage hmor
  refine ⟨f.left.hom, f.left.property.1, f.left.property.2, ?_⟩
  have hw := Over.w f
  exact congrArg (fun t => t.hom) hw

include P in
/-- ★★**第1の圏同値(コスライス側)を「因子分解」として消費する**。

`coaPre_factor_of_mle` の**鏡像**である —— あちらは `Cc` へ**入る** co-angular
pre-step どうしを比べ、こちらは `A` から**出る**ものどうしを比べる。

★**向きの違いが 2 か所に出る**:
- 不変量が `Φ.map (inv (Base γ)) (Div γ)` ではなく **`Div ψ` そのもの**
  (コスライスでは `Φ(A)` が既に始域側の単系だから、輸送が要らない)
- `Order` への関手が**共変**なので、`op` を取らない

★★**`Proposition 1.14, (i)` の「`Div(φ)` が irreducible」がこれを要求する** ——
`Div φ = a + b` という分解から `φ` の分解を作るのに、
「`a` を実現する `ψ`(`coaPre_realize`)＋ `φ` が `ψ` を経由すること(この補題)」の
2 つが要る。 -/
theorem coaPre_factor_under_of_mle (G : Frobenioid P) {A X B : C}
    (ψ : A ⟶ X) (hψc : IsCoAngular P ψ) (hψs : IsPreStep P ψ)
    (φ : A ⟶ B) (hφc : IsCoAngular P φ) (hφs : IsPreStep P φ)
    (hle : MLe (P.Div ψ) (P.Div φ)) :
    ∃ β : X ⟶ B, IsCoAngular P β ∧ IsPreStep P β ∧ ψ ≫ β = φ := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  set Zψ : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨X⟩ from ⟨ψ, hψc, hψs⟩) with hZψ
  set Zφ : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨B⟩ from ⟨φ, hφc, hφs⟩) with hZφ
  have hmor : (coaPreUnderFunctor P A).obj Zψ ⟶ (coaPreUnderFunctor P A).obj Zφ :=
    homOfLE (show toOrderCat (P.Div ψ) ≤ toOrderCat (P.Div φ) from hle)
  let f := (coaPreUnderFunctor P A).preimage hmor
  refine ⟨f.right.hom, f.right.property.1, f.right.property.2, ?_⟩
  exact congrArg (fun t => t.hom) (Under.w f)

/-! ### ★(vi) の前半の組み立て

原文(p.36)の段取りは 7 段:
1. `Definition 1.3, (i), (a), (b)` で co-angular pre-step `α : B ⟶ A`、`γ : B ⟶ C`
   (`C` は Frobenius-trivial)を得る
2. `C` が Frobenius-trivial なので、各 `d` に base-identity 自己射 `φ_C`(次数 `d`)がある
3. **(ii)** で `γ ≫ φ_C = ψ ≫ γ′` と組み替える(`ψ` は Frobenius 型、`γ′` は co-angular pre-step)
4. **`Div` の関係 ＋ 圏同値**で `γ′` が `γ` を経由する: `β ≫ γ = γ′`
5. `φ_B := ψ ≫ β` と置くと `φ_B ≫ γ = γ ≫ φ_C`
6. `prop_1_10_vi_descend` で `φ_B` は base-identity・次数 `d`
7. よって `B` は quasi-Frobenius-trivial、`α` により `A` は sub-quasi-Frobenius-trivial

★★**ここでは 6・7 を実装し、3–5 が作るべきものを `hstep` として仮定に出す。**
★**こうすると「残りは何か」が宣言の型として確定する** ——
`Prop110iGoal` で使ったのと同じやり方である。
★**未完を散文で書くのではなく、型で書く。**
-/

include P in
/-- **(vi) 前半の 6・7 段** —— `γ` に沿って base-identity 自己射が降りるなら
`B` は quasi-Frobenius-trivial。

★仮定 `hstep` が原文の 3–5 段(**(ii) と圏同値**)が作るべきものである。
★**`prop_1_10_vi_descend` がそれを受けて 6 段を実行する。** -/
theorem prop_1_10_vi_quasi {B Cc : C} (γ : B ⟶ Cc) (hγs : IsPreStep P γ)
    (hC : IsFrobeniusTrivial P Cc)
    (hstep : ∀ (n : ℕ+) (φC : Cc ⟶ Cc), IsBaseIdentity P φC → IsFrobeniusType P φC →
      P.degFr φC = n → ∃ φB : B ⟶ B, φB ≫ γ = γ ≫ φC) :
    IsQuasiFrobeniusTrivial P B := by
  intro n
  obtain ⟨ζ, hdeg, hprop⟩ := hC
  obtain ⟨φB, hsq⟩ := hstep n (ζ n) (hprop n).1 (hprop n).2 (hdeg n)
  obtain ⟨hb, hd⟩ := prop_1_10_vi_descend P γ hγs (ζ n) φB (hprop n).1 (hdeg n) hsq
  exact ⟨φB, hb, hd⟩

include P in
/-- **(vi) 前半の 7 段** —— quasi-Frobenius-trivial な対象へ co-angular pre-step が
あれば sub-quasi-Frobenius-trivial。

★**定義そのものだが、原文の「hence that A is sub-quasi-Frobenius-trivial」の
一言に対応する**ので明示する。 -/
theorem prop_1_10_vi_subQuasi {A B : C} (α : B ⟶ A)
    (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (hq : IsQuasiFrobeniusTrivial P B) : IsSubQuasiFrobeniusTrivial P A :=
  ⟨B, α, hαc, hαs, hq⟩

/-! ### ★MLe の補助 —— 穴を塞ぐのに要る

★`MLe a b := ∃ c, a + c = b` の定義から直に出る。
★**`mle_nsmul_self` は 2026-08-17 に `MonoidPrime.lean` へ移した**——
§0 の語彙なのだから、実装場所もそちらであるべきだった。
-/

/-- ★**加法準同型は `MLe` を保つ**。 -/
theorem MLe.map {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    {a b : M} (h : MLe a b) : MLe (f a) (f b) := by
  obtain ⟨c, hc⟩ := h
  exact ⟨f c, by rw [← map_add, hc]⟩

/-- ★**加法準同型は `≼` も保つ** —— `MLe.map` に `map_nsmul` を挟むだけ。

★★**`Proposition 1.14, (v)` で `non-dilating` の仮定を作るのに要る。** -/
theorem MPrec.map {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    {a b : M} (h : MPrec a b) : MPrec (f a) (f b) := by
  obtain ⟨n, hn, hle⟩ := h
  exact ⟨n, hn, by simpa [map_nsmul] using MLe.map f hle⟩

/-! ### ★★(vi) の最後の穴を塞ぐ —— 3–5 段

★**道具は全部そろっている**: (ii)(`prop_1_10_ii`)、`Div` の公式
(`prop_1_10_ii_Div_formula`)、橋(`coaPre_factor_of_mle`)。
★**残っていたのは `Φ.map` の記帳だけ**だった。

★**記帳の中身**:
- 四角形 `β′ ≫ α′ = γ ≫ φ_C` から `Base γ = Base β′ ≫ Base α′`(`φ_C` は base-identity)
- `Div` の公式から `Φ.map (Base β′) (Div α′) = d • Div γ`
- 両辺に `Φ.map (inv (Base β′))` を当てて ★**`Div α′ = d • Φ.map (inv (Base β′)) (Div γ)`**
- `x ≤ d • x` なので `MLe (Φ.map (inv (Base β′)) (Div γ)) (Div α′)`
- `Φ.map (inv (Base α′))` で押し出すと、★**橋が要求する形**になる
-/

include P in
/-- ★★**(vi) 前半の 3–5 段** —— `prop_1_10_vi_quasi` が要求する `hstep` を作る。

★これで `Proposition 1.10, (vi)` の前半の穴が塞がる。 -/
theorem prop_1_10_vi_step (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {B Cc : C} (γ : B ⟶ Cc) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (n : ℕ+) (φC : Cc ⟶ Cc) (hbC : IsBaseIdentity P φC) (hftC : IsFrobeniusType P φC)
    (hdC : P.degFr φC = n) :
    ∃ φB : B ⟶ B, φB ≫ γ = γ ≫ φC := by
  -- 3 段: (ii) で組み替える
  obtain ⟨Y, β', α', hβ', hdβ', hα's, hsq⟩ := prop_1_10_ii P G.core γ hγs φC hftC
  haveI hbβ' : IsIso (P.Base β') := hβ'.2
  haveI hbα' : IsIso (P.Base α') := hα's.2
  haveI hbγ : IsIso (P.Base γ) := hγs.2
  -- `Div` の公式
  have hdiv : Φ.map (P.Base β') (P.Div α') = (P.degFr φC : ℕ) • P.Div γ :=
    prop_1_10_ii_Div_formula P γ φC β' α' hβ' hftC hsq
  -- `Div α′ = d • Φ.map (inv (Base β′)) (Div γ)`
  have hdiv' : P.Div α' = (P.degFr φC : ℕ) • Φ.map (inv (P.Base β')) (P.Div γ) := by
    have := congrArg (Φ.map (inv (P.Base β'))) hdiv
    rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id, map_nsmul] at this
    exact this
  -- `MLe` を作る
  have hmle0 : MLe (Φ.map (inv (P.Base β')) (P.Div γ)) (P.Div α') := by
    rw [hdiv']
    exact mle_nsmul_self (P.degFr φC).pos _
  -- 橋が要求する形へ押し出す
  -- ★`inv` の下を書き換えると motive エラー(分類表 #2)になるので、
  --   `inv_eq_of_hom_inv_id` で**外から**特徴づける
  have hb : P.Base β' ≫ P.Base α' = P.Base γ := by
    have h := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp, show P.Base φC = P.Base (𝟙 Cc) from hbC,
      P.Base_id, Category.comp_id] at h
    exact h
  have hbase : inv (P.Base γ) = inv (P.Base α') ≫ inv (P.Base β') := by
    refine CategoryTheory.IsIso.inv_eq_of_hom_inv_id ?_
    rw [← hb, Category.assoc, ← Category.assoc (P.Base α'), IsIso.hom_inv_id,
      Category.id_comp, IsIso.hom_inv_id]
  have hmle : MLe (Φ.map (inv (P.Base γ)) (P.Div γ))
      (Φ.map (inv (P.Base α')) (P.Div α')) := by
    rw [hbase, Φ.map_comp]
    exact MLe.map _ hmle0
  -- 4 段: 橋で因子分解を得る
  obtain ⟨β, _, _, hβγ⟩ :=
    coaPre_factor_of_mle P G γ hγc hγs α' (prop_1_4_i P α' (fun Y' _ => hiso Y')) hα's hmle
  -- 5 段: `φ_B := β′ ≫ β`
  exact ⟨β' ≫ β, by rw [Category.assoc, hβγ, hsq]⟩

include P in
/-- ★★**`Proposition 1.10, (vi)` 前半の完成形** ——
isotropic 型の Frobenioid において、Frobenius-trivial な対象へ co-angular pre-step で
つながる対象は quasi-Frobenius-trivial であり、そこへ co-angular pre-step で
つながる対象は sub-quasi-Frobenius-trivial。

★原文(p.36)の 7 段をすべて実装した形である。
★**`Definition 1.3, (i), (a), (b)` から `α`・`γ`・`C` を得る 1・2 段は仮定に出している**
——それは `𝒞^istr` が Frobenioid であること(`Proposition 1.9, (v)`)から供給されるもので、
**この定理の内容ではない**。 -/
theorem prop_1_10_vi_first (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B Cc : C} (α : B ⟶ A) (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (γ : B ⟶ Cc) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (hC : IsFrobeniusTrivial P Cc) :
    IsSubQuasiFrobeniusTrivial P A :=
  prop_1_10_vi_subQuasi P α hαc hαs
    (prop_1_10_vi_quasi P γ hγs hC
      (fun n φC hbC hftC hdC =>
        prop_1_10_vi_step P G hiso γ hγc hγs n φC hbC hftC hdC))

include P in
/-- ★★★**`α`・`γ`・`C` を `Definition 1.3, (i), (a), (b)` から導く**（2026-08-16）。

★★**監査で「原文は導いているのに仮定に逃がしている」と指摘された段**である。

★**導き方**（3 段）:
1. `baseSurj` を `Base A` に当て、**Frobenius-trivial な `A₀` と底の同型**を得る
2. `preStepSpan` をその同型に当て、**pre-step の span** `A ← X → A₀` を得る
3. ★`𝒞` が isotropic 型なので、`Proposition 1.4, (i)` により
   **すべての射が co-angular** —— span の両辺が co-angular pre-step になる

★★**これで `prop_1_10_vi_first` の仮定がすべて供給される**。 -/
theorem prop_1_10_vi_data (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X) (A : C) :
    ∃ (B Cc : C) (α : B ⟶ A) (γ : B ⟶ Cc),
      IsCoAngular P α ∧ IsPreStep P α ∧ IsCoAngular P γ ∧ IsPreStep P γ ∧
      IsFrobeniusTrivial P Cc := by
  obtain ⟨A₀, hA₀ft, ⟨e⟩⟩ := G.core.baseSurj (P.toElem.obj A).base
  haveI : IsIso e.inv := inferInstance
  obtain ⟨X, φ, ψ, hφ, hψ, -⟩ := G.core.preStepSpan A₀ A e.hom e.isIso_hom
  exact ⟨X, A₀, ψ, φ,
    prop_1_4_i P ψ (fun Y _ => hiso Y), hψ,
    prop_1_4_i P φ (fun Y _ => hiso Y), hφ, hA₀ft⟩

include P in
/-- ★★★**`Proposition 1.10, (vi)` 前半の「型」の形**（2026-08-16）。

原文 (FrdI p.35):
> (vi) The Frobenioid Cistr is of sub-quasi-Frobenius-trivial type. Moreover,

★★**監査で「「type」の ∀ が無い」と指摘された形**である。
`prop_1_10_vi_data`（`Definition 1.3, (i), (a), (b)` からの導出）と
`prop_1_10_vi_first`（本体）を合成するだけで出る。

★**isotropic 型の Frobenioid は sub-quasi-Frobenius-trivial 型である** ——
★`𝒞^istr` は isotropic 型なので、これを `istrPre` に当てれば原文の形になる。 -/
theorem prop_1_10_vi_ofType (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X) :
    IsOfSubQuasiFrobeniusTrivialType P := by
  intro A
  obtain ⟨B, Cc, α, γ, hαc, hαs, hγc, hγs, hC⟩ := prop_1_10_vi_data P G hiso A
  exact prop_1_10_vi_first P G hiso α hαc hαs γ hγc hγs hC


/-! ## ★★(iii) —— 第1の圏同値も消費する

★**前段では第2(スライス)の圏同値を消費した**(`coaPre_factor_of_mle`)。
(iii) では **第1(コスライス)の圏同値**を使う。
関手は `_A(𝒞^coa-pre) ⥤ Order(Φ(A))`、`Z ↦ Div (Z.hom)`。

★**本質的全射性**は「任意の `c ∈ Φ(A)` は、ある co-angular pre-step の `Div` に**同型**」
を与える。★**`Order` は前順序圏なので、同型 = 両向きの `MLe`。**
★**そこに `mle_antisymm`(integral ＋ sharp)を当てると等号になる。**
★★**`divisorial` の 2 つの条件(integral と sharp)が、ここで一緒に効く。**
-/

include P in
/-- ★★**第1の圏同値の初の消費** —— `Φ(A)` の任意の元は
co-angular pre-step の `Div` として実現できる。

★**前順序圏の同型は両向きの `MLe`** なので、`mle_antisymm` で等号に直す。 -/
theorem coaPre_realize (G : Frobenioid P) (A : C)
    (c : Φ.val (P.toElem.obj A).base) :
    ∃ (X : C) (ψ : A ⟶ X), IsCoAngular P ψ ∧ IsPreStep P ψ ∧ P.Div ψ = c := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := coaPreUnderFunctor P A)
    (toOrderCat c)
  refine ⟨Z.right.obj, Z.hom.hom, Z.hom.property.1, Z.hom.property.2, ?_⟩
  refine mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 ?_ ?_
  · exact leOfHom e.hom
  · exact leOfHom e.inv

/-! ### ★(iii) の全射性 —— `n` 倍が全射であること

★**組み立て**:
1. `c ∈ Φ(A)` を `coaPre_realize` で co-angular pre-step `ψ : A ⟶ X`(`Div ψ = c`)に実現
2. `IsPerfectObj A` の (a) で、`A` と `X` の上に**次数 `n` の Frobenius 型射**
   `φ₁ : B₀ ⟶ A`、`φ₂ : B₂ ⟶ X` を取る(`ψ` が pre-step なので `X` は `A` と base-isomorphic)
3. `IsPerfectObj A` の (b) で `ψ` を降ろし、pre-step `ψ₀ : B₀ ⟶ B₂` と
   四角形 `φ₁ ≫ ψ = ψ₀ ≫ φ₂` を得る
4. ★**核心の等式**(`prop_1_10_iii_div_nsmul'`)から `Φ.map (Base φ₁) (Div ψ) = n • Div ψ₀`
5. `Φ.map (inv (Base φ₁))` を当てて ★**`c = n • Φ.map (inv (Base φ₁)) (Div ψ₀)`**

★**「A は perfect」から「Φ(A) は perfect」への橋が、これで渡り切る。**
-/

include P in
/-- ★★**(iii) の全射性** —— `𝒞` の対象が perfect なら、`Φ(A)` で `n` 倍は全射。 -/
theorem prop_1_10_iii_nsmul_surjective (G : Frobenioid P) {A : C}
    (hperf : IsPerfectObj P A) (n : ℕ+) :
    Function.Surjective (fun a : Φ.val (P.toElem.obj A).base => (n : ℕ) • a) := by
  intro c
  obtain ⟨X, ψ, hψc, hψs, hψd⟩ := coaPre_realize P G A c
  haveI hbψ : IsIso (P.Base ψ) := hψs.2
  have hAX : BaseIsomorphic P A X := ⟨asIso (P.Base ψ)⟩
  have hAA : BaseIsomorphic P A A := ⟨Iso.refl _⟩
  obtain ⟨B₀, φ₁, hφ₁, hd₁⟩ := (hperf n).1 A hAA
  obtain ⟨B₂, φ₂, hφ₂, hd₂⟩ := (hperf n).1 X hAX
  haveI hbφ₁ : IsIso (P.Base φ₁) := hφ₁.2
  -- `B₀`・`B₂` は `A` と base-isomorphic
  have hB₀A : BaseIsomorphic P B₀ A := ⟨asIso (P.Base φ₁)⟩
  have hB₂A : BaseIsomorphic P B₂ A := by
    haveI : IsIso (P.Base φ₂) := hφ₂.2
    exact ⟨asIso (P.Base φ₂) ≪≫ (asIso (P.Base ψ)).symm⟩
  obtain ⟨ψ₀, ⟨hψ₀s, hsq⟩, _⟩ :=
    (hperf n).2 B₀ A B₂ X φ₁ φ₂ hφ₁ hd₁ hφ₂ hd₂ hB₀A hB₂A ψ hψs
  -- 核心の等式
  have hcore : Φ.map (P.Base φ₁) (P.Div ψ) = (n : ℕ) • P.Div ψ₀ :=
    prop_1_10_iii_div_nsmul' P φ₁ φ₂ ψ₀ ψ hφ₁.1.2 hφ₂.1.2 hd₂ hsq
  refine ⟨Φ.map (inv (P.Base φ₁)) (P.Div ψ₀), ?_⟩
  have h := congrArg (Φ.map (inv (P.Base φ₁))) hcore
  rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id, map_nsmul] at h
  show (n : ℕ) • Φ.map (inv (P.Base φ₁)) (P.Div ψ₀) = c
  rw [← h, hψd]

/-! ### ★★(iii) の単射性 —— **親の前段の見立ては誤っていた。そして次の壁が見えた**

前段で親はこう書いた:
> 単射性は `IsPerfectObj` の `∃!` から出るはずで、`divisorial` からは出ない
> (`Gp M` の捻れを排除する条件が無い)。

★★**誤りである。紙の上では `divisorial` から出る。** `Gp M` は**捻れなし**である:

`x ∈ Gp M` が `n • x = 0`(`n ≥ 1`)を満たすとする。
- `n • x = 0 = toGp 0` は `M` の像に入るので、★**saturated** から `x = toGp m` と書ける
- 同様に `n • (-x) = 0` も像に入るので `-x = toGp m′`
- `toGp (m + m′) = x + (-x) = 0 = toGp 0` で、★**integral**(`toGp` が単射)から `m + m′ = 0`
- よって `m` は可逆、★**sharp** から `m = 0`、したがって `x = 0`

★★**`divisorial` の 3 条件(saturated / integral / sharp)が、この 1 つの主張で
一度に効いている。**

## ★★しかし Lean では書けなかった —— 次の壁

★**`Gp M := AddLocalization (⊤ : AddSubmonoid M)` に群構造が無い。**
上の議論は `x - y` と `-x` を使うが、我々の `Gp M` は
**`AddCommMonoid` としてしか実装されていない**(`HSub`・`Neg` のインスタンスが無い)。
★**`saturated` の定義自体が `Gp M` の元を量化しているのに、
`Gp M` を群として使う準備ができていない。**

★**したがって単射性はここで止まる。** `sorry` は置かない(`Found/` の規律)。
★**次に当たるときの仕事は明確**: `Gp M` に `AddCommGroup` インスタンスを付ける
(`neg (mk a s) := mk s a`、`S = ⊤` なので常に定義できる)。

## ★★記録: 否定的な予測は検証が甘くなる

親は前段で「`divisorial` からは出ない」と書いた。★**1 本も試さずに書いた。**
実際に試したら**紙の上では出た**(そして Lean では別の理由で止まった)。

★★**否定的な予測(「〜からは出ない」)は、肯定的な予測より検証が甘くなる。**
「出る」は 1 本示せば済むが、「出ない」は**全ての道を塞いだこと**を要する。
★**我々は「証明できなかった」と「反例がある」を区別してきたが、
「試していない」と「試して駄目だった」の区別も要る。**
★今回の親の発言は **3 つ目のカテゴリ(試していない)**だった。
-/

/-! ### ★★壁は import 1 行だった —— `Gp M` の群構造

前段で親は「`Gp M` に群構造が無いので単射性は書けない」と記録した。
★**確かめたら、mathlib に**あった**。**

`Mathlib/GroupTheory/MonoidLocalization/GrothendieckGroup.lean`:
```
abbrev AddGrothendieckGroup : Type _ := AddLocalization (⊤ : AddSubmonoid M)
instance instAddCommGroup : AddCommGroup (AddGrothendieckGroup M)
```
★★**我々の `Gp M := AddLocalization (⊤ : AddSubmonoid M)` と同じもの**である。
**import していなかっただけ。**

★**これは `comp_div` と同じ形の壁である**(第24段): 症状は「難しい」に見えたが、
原因は★**「1 行足りない」**だった。
★**前回は `@[simp]` 補題が 1 本、今回は `import` が 1 行。**
★★**「壁の高さ」と「原因の大きさ」は関係がない。**

★**そして 2 回とも、原因は「揃っているはずのものを並べて見る」で見つかった** ——
前回は `@[simp]` 補題の一覧、今回は mathlib のファイル一覧
(`ls` して `GrothendieckGroup.lean` が目に入った)。
-/

/-- ★★**divisorial なモノイドでは `n` 倍が単射**(`n ≥ 1`)。

★`Gp M` が**捻れなし**であることの帰結。saturated / integral / sharp の 3 つが要る。 -/
theorem nsmul_injective_of_isDivisorial {M : Type w} [AddCommMonoid M]
    (hd : IsDivisorial M) {n : ℕ} (hn : 0 < n) :
    Function.Injective (fun a : M => n • a) := by
  obtain ⟨⟨hint, hsat, _⟩, hsh⟩ := hd
  -- ★`AddLocalization.mk_add`(mathlib、`mk_mul` の `to_additive` 版)で足し算を開く
  have htoGp_add : ∀ x y : M, toGp M (x + y) = toGp M x + toGp M y := by
    intro x y
    rw [toGp, toGp, toGp, AddLocalization.mk_add]
    congr 1
    exact (add_zero (0 : (⊤ : AddSubmonoid M))).symm
  have htoGp_nsmul : ∀ (k : ℕ) (x : M), toGp M (k • x) = k • toGp M x := by
    intro k x
    induction k with
    | zero => simp [toGp, AddLocalization.mk_zero]
    | succ m ih => rw [succ_nsmul, succ_nsmul, htoGp_add, ih]
  intro a b hab
  set x : Gp M := toGp M a - toGp M b with hx
  have hnx : n • x = 0 := by
    rw [hx, smul_sub, ← htoGp_nsmul, ← htoGp_nsmul, sub_eq_zero]
    exact congrArg (toGp M) hab
  have hmem : ∀ y : Gp M, n • y = 0 → y ∈ Set.range (toGp M) := by
    intro y hy
    refine hsat y n hn ?_
    rw [hy]
    exact ⟨0, by simp [toGp, AddLocalization.mk_zero]⟩
  obtain ⟨m, hm⟩ := hmem x hnx
  obtain ⟨m', hm'⟩ := hmem (-x) (by rw [smul_neg, hnx, neg_zero])
  have hsum : m + m' = 0 := by
    refine hint ?_
    rw [htoGp_add, hm, hm', add_neg_cancel]
    simp [toGp, AddLocalization.mk_zero]
  have hm0 : m = 0 := hsh m ⟨⟨m, m', hsum, by rw [add_comm]; exact hsum⟩, rfl⟩
  have hx0 : x = 0 := by rw [← hm, hm0]; simp [toGp, AddLocalization.mk_zero]
  refine hint ?_
  rwa [hx, sub_eq_zero] at hx0

include P in
/-- ★★**`Proposition 1.10, (iii)` の第1主張 完成** ——
`𝒞` の対象 `A` が perfect なら、`Φ(A)` は perfect なモノイド。

★**全射性**は `IsPerfectObj`(原文が名指す `Definition 1.2, (iv)`)から、
★**単射性**は `divisorial`(saturated ＋ integral ＋ sharp)から。

★★**原文は「the fact that Φ(A) is perfect follows immediately … from the fact that
A is perfect」と書くが、実際には 2 つの別々の源から来ている** ——
全射性だけが `A` の perfect 性から来て、**単射性は `Φ` が divisorial であることから来る**。
★**原文の「from the fact that A is perfect」は半分しか説明していない。** -/
theorem prop_1_10_iii_perfect (G : Frobenioid P) {A : C} (hperf : IsPerfectObj P A) :
    IsPerfectMonoid (Φ.val (P.toElem.obj A).base) := fun n =>
  ⟨nsmul_injective_of_isDivisorial (P.divisorial _) n.pos,
   prop_1_10_iii_nsmul_surjective P G hperf n⟩

include P in
/-- ★**(iii) 第1主張の「型」版** —— `𝒞` が of perfect type なら
`Φ` の像のモノイドはすべて perfect。★**原文の言い回しそのもの。** -/
theorem prop_1_10_iii_image_perfect (G : Frobenioid P) (hpt : IsOfPerfectType P) (A : C) :
    IsPerfectMonoid (Φ.val (P.toElem.obj A).base) :=
  prop_1_10_iii_perfect P G (hpt A)

include P in
/-- ★★**(iii) 第1主張 ——「the monoids in the image of Φ」の本当の範囲**。

★★**取り下げの原因の1つを埋める**(2026-08-16)。
`prop_1_10_iii_image_perfect` は `𝒞` の対象の底、すなわち
**`Base` の像**についてしか述べていなかった。★原文の「the image of Φ」は
`Φ : 𝒟ᵒᵖ ⥤ 𝔐𝔬𝔫` の像であり、★**`Ob(𝒟)` 全体に渡る。**

★**渡し方**: `Definition 1.3, (i), (a)`(`baseSurj`)が、任意の `Y : 𝒟` に対し
`Base(A) ≅ Y` なる(しかも Frobenius-trivial な)`A : 𝒞` を与える。
`Φ.map` はその同型に沿って全単射(`MonoidOn.map_bijective_of_iso`)であり、
perfect 性は加法的全単射で移る(`isPerfectMonoid_of_bijective`)。

★★**ここで初めて `𝒟` 側と `𝒞` 側の隔たりが埋まった** ——
`Proposition 1.10` の他の条はすべて `𝒞` の中の話だが、
★**(iii) の第1主張だけは `𝒟` の全対象について述べている。**
その差を原文は書いていない(「the image of Φ」の一語で済ませている)。 -/
theorem prop_1_10_iii_image_perfect_of_base (F : FrobenioidCore P) (G : Frobenioid P)
    (hpt : IsOfPerfectType P) (Y : D) : IsPerfectMonoid (Φ.val Y) := by
  obtain ⟨A, -, ⟨e⟩⟩ := F.baseSurj Y
  exact isPerfectMonoid_of_bijective (Φ.map e.symm.hom) (Φ.map_bijective_of_iso e.symm)
    (prop_1_10_iii_image_perfect P G hpt A)

end ABC3.Found.FrdI
