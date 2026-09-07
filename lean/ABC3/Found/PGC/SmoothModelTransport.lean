import ABC3.Found.PGC.DegreeTransport
import ABC3.Skeleton.PGC.Section2Defs
import Mathlib.Topology.LocallyConstant.Algebra
import Mathlib.FieldTheory.Galois.NormalBasis

/-!
# `Γ_K`-加群の群論的移送 —— [pGC] Proposition 2.1 への抽象核と還元

原文 (pGC 物理 p.4, Proposition 2.1):

> The Γ_K-module K[bar] may be recovered group-theoretically from Γ_K.

本ファイルは `Skeleton/PGC/Section2.lean::prop_2_1` を

  「`K̄` は `Γ_K` の位相群としての構造だけから作れる加群と同変同型である」

という**ただ 1 つの仮説**に還元する。還元そのものは `sorry` 無しである。

## 1. 抽象核(★分岐・付値・Galois の語彙が 1 語も出ない)

* `SemilinearAddEquiv α A B` —— 群同型 `α : G ≃* G'` に沿って同変な加法同型。
  `refl`/`symm`/`trans`/`ofEq`/`pi`/`sandwich` を備える。純群論+加法群のみ。
* `locallyConstantTranslation` —— 位相群 `G` 上の局所定数関数 `LocallyConstant G M`
  への平行移動作用 `(g • f) x = f (g⁻¹ x)`。
* `locallyConstantSourceCongr` —— `α : G ≃ₜ* G'` から
  `SemilinearAddEquiv α.toMulEquiv (LocallyConstant G M) (LocallyConstant G' M)`。
  ★**これが「群論的に回復できる」の内実**である:`LocallyConstant G M` は
  `G` の位相群としての構造と `M` だけから作られるので、`α` がそのまま同型を運ぶ。

## 2. 具体層 —— 残った 1 つの仮説

`SmoothModel K d` を

  `K̄ ≅ LocallyConstant Γ_K (ℚ_p^d)`  (`Γ_K`-同変)

と置くと、`prop_2_1_of_smoothModel` が

  `(∀ K, SmoothModel K [K:ℚ_p]) → RecoverableAsAddModule (fun K => K.closure)`

を与える。`[K:ℚ_p] = [K':ℚ_p]` は在庫
`Found/PGC/DegreeTransport.lean::finrank_eq_of_absGal_equiv`(`sorry` 無し)から出る。

`SmoothModel` は**古典的な「無限次 Galois 拡大に対する正規底定理」**そのものである:
`K̄ = colim_L L`、有限次 Galois `L/K` に対し正規底定理が `L ≅ K[Gal(L/K)]`
を与え、その両立系の colimit が `LocallyConstant Γ_K K` になる
(基底元の移行は跡写像 `Tr_{L'/L}` で、`K[G']^H ≅ K[G]` の下で基底元は基底元に落ちる)。
`K ≅ ℚ_p^{[K:ℚ_p]}` を入れれば上の形になる —— それが `smoothModel_of_carrier`。

★★**有限次の正規底定理は mathlib に在る**(実測 2026-09-08、
`grep -n "NormalBasis" .cache/mathlib-index.txt | grep -v -i orthonormal`):

    Mathlib/FieldTheory/Galois/NormalBasis.lean
      IsGalois.normalBasis (K L) : Module.Basis Gal(L/K) K L
      IsGalois.normalBasis_apply (e : Gal(L/K)) : normalBasis K L e = e (normalBasis K L 1)

★**私は最初これを「無い」と誤判定した**——`grep -i "normalBasis"` は
`OrthonormalBasis` を大量に拾い、`head -30` で切れて `FieldTheory/Galois/NormalBasis.lean`
が視野の外に出ていた。`lean-idioms.md` #117(ii) の失敗形そのものである。

本ファイル §5 はその mathlib の基底を**同変形に上げる**
(`finiteNormalBasisModel`:有限次 Galois `E/F` に対し
`E ≅ LocallyConstant Gal(E/F) F` が `Gal(E/F)`-同変、`sorry` 無し)。
mathlib は「基底がある」までしか言わないので、同変性は本ファイルの新しい内容である。

★残った穴は**有限次から無限次への coherent な極限**のみ:
`L ⊆ L'` に対し正規底生成元 `θ'` の跡 `Tr_{L'/L} θ'` は `L` の正規底生成元になる
(`K[G']^H ≅ K[G]` の下で基底元は基底元に落ちる)ので降下は自動、
上げる側は `K[G'] ↠ K[G]` の単元の持ち上げ(Maschke で `K[G']` が半単純なので分裂)で済む。
そこを formalize する節点は本ファイルには含まれていない。

## ★★原典より短い道(記録)

原典 §2 の論拠は p 進対数 + Verlagerung + 局所類体論である:

> Recall that the p-adic logarithm defines a natural isomorphism of UK (modulo torsion)
> onto an open subgroup of K. In particular, it defines an isomorphism of UK ⊗Zp Qp with
> K. … if we regard UK (respectively, UL) as a subgroup of Γab_K (respectively, Γab_L)
> then the morphism UK →UL may be recovered group-theoretically by means of the
> "Verlagerung, or transfer, map"

本ファイルが採るのは別の道であり、次の 3 つを一切使わない:

1. p 進対数(在庫はある —— `Found/PGC/PrincipalUnitsLog.lean::padicLogUnitsEquiv`、
   `Found/PGC/PrincipalUnitsRank.lean::finrank_smallBall`)。
2. Verlagerung(mathlib `Mathlib/GroupTheory/Transfer.lean::MonoidHom.transfer`)。
3. 局所類体論の相互写像(`Found/PGC/ArtinMap.lean::artinMapAbelian`)。

★理由(実測):原典の道は `U_K ⊆ Γ_K^ab` の**同変性**(Artin 写像が `α` と可換であること)
を要求する。それは本木では `Found/PGC/ArtinEquivariance.lean` として**未解決のまま
残っている**(`ArtinEquivariance`、`sorry` 3 件)。さらに「`U_K → U_L` が Verlagerung と
一致する」([6] Serre VII §8)は mathlib にも本木にも無い。
一方こちらの道が要求するのは正規底定理 1 本だけで、`[K:ℚ_p]` の移送は既に
`sorry` 無しで在庫にある。

★★これは逸脱ではない:得られる `RecoverableAsAddModule` の**主張は同一**であり、
証明の道筋だけが違う(原典が `U_K ⊗ Q_p ≅ K` を経由するところを、
正規底 `K̄ ≅ C^∞(Γ_K, K)` で置き換えた)。

## 3. 非空虚性(★仮説を空虚に満たしていないことの検査)

`recoverableAsAddModule_locallyConstant` は、`Obj K := LocallyConstant Γ_K M`
という族について `RecoverableAsAddModule` が**無条件に成り立つ**ことを示す
(`sorry` 無し)。すなわち §2 の還元先は空虚ではない。
また `smul_eq_self_iff_const` は `(LocallyConstant G M)` の `G`-不変元が
ちょうど定数関数であることを示す —— `K̄^{Γ_K} = K` と合わせると、模型の係数は
`ℚ_p^{[K:ℚ_p]}` 以外にありえない。すなわち `SmoothModel` の `d` は
勝手に選べる自由度ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

universe u v w

/-! ## 1. ★★★抽象核 I —— `α` に沿って同変な加法同型(語彙ゼロ) -/

/-- 群同型 `α : G ≃* G'` に沿って**同変**な加法群の同型。

`RecoverableAsAddModule` の結論部
`∃ φ : Obj K ≃+ Obj K', ∀ g x, φ (g • x) = α g • φ x` を、
局所体を一切含まない形で切り出したもの。 -/
structure SemilinearAddEquiv {G : Type u} {G' : Type v} [Group G] [Group G']
    (α : G ≃* G') (A : Type w) (B : Type*)
    [AddCommGroup A] [AddCommGroup B] [DistribMulAction G A] [DistribMulAction G' B] where
  /-- 台となる加法同型 -/
  toAddEquiv : A ≃+ B
  /-- `α` に沿った同変性 -/
  map_smul : ∀ (g : G) (a : A), toAddEquiv (g • a) = α g • toAddEquiv a

namespace SemilinearAddEquiv

variable {G G' G'' : Type*} [Group G] [Group G'] [Group G'']
variable {A B C D : Type*} [AddCommGroup A] [AddCommGroup B] [AddCommGroup C] [AddCommGroup D]
variable [DistribMulAction G A] [DistribMulAction G' B] [DistribMulAction G'' C]

/-- 恒等。 -/
def refl : SemilinearAddEquiv (MulEquiv.refl G) A A :=
  ⟨AddEquiv.refl A, fun _ _ => rfl⟩

/-- 逆。 -/
def symm {α : G ≃* G'} (e : SemilinearAddEquiv α A B) : SemilinearAddEquiv α.symm B A :=
  ⟨e.toAddEquiv.symm, by
    intro g' b
    apply e.toAddEquiv.injective
    rw [e.map_smul]
    simp⟩

/-- 合成。 -/
def trans {α : G ≃* G'} {β : G' ≃* G''} (e : SemilinearAddEquiv α A B)
    (f : SemilinearAddEquiv β B C) : SemilinearAddEquiv (α.trans β) A C :=
  ⟨e.toAddEquiv.trans f.toAddEquiv, by
    intro g a
    simp only [AddEquiv.trans_apply, MulEquiv.coe_trans, Function.comp_apply]
    rw [e.map_smul, f.map_smul]⟩

/-- `α` を等しいものに取り替える。 -/
def ofEq {α α' : G ≃* G'} (h : α = α') (e : SemilinearAddEquiv α A B) :
    SemilinearAddEquiv α' A B := by
  subst h; exact e

variable {ι : Type*}

/-- 添字づけた直積へ持ち上げる。 -/
def pi {α : G ≃* G'} (e : SemilinearAddEquiv α A B) :
    SemilinearAddEquiv α (ι → A) (ι → B) :=
  ⟨AddEquiv.piCongrRight (fun _ : ι => e.toAddEquiv), by
    intro g a
    funext i
    exact e.map_smul g (a i)⟩

/-- ★**「模型を経由する」を 1 本にまとめた核**。

`A` が `G`-同変に `C` と、`B` が `G'`-同変に `D` と同型で、
`C` と `D` が `α` に沿って同型なら、`A` と `B` が `α` に沿って同型。
[pGC] Proposition 2.1 の証明はこの形をしている(`C`・`D` が「群論的な模型」)。 -/
def sandwich {α : G ≃* G'} [DistribMulAction G C] [DistribMulAction G' D]
    (eA : SemilinearAddEquiv (MulEquiv.refl G) A C)
    (e : SemilinearAddEquiv α C D)
    (eB : SemilinearAddEquiv (MulEquiv.refl G') B D) :
    SemilinearAddEquiv α A B :=
  ⟨eA.toAddEquiv.trans (e.toAddEquiv.trans eB.toAddEquiv.symm), by
    intro g a
    have h1 : eA.toAddEquiv (g • a) = g • eA.toAddEquiv a := eA.map_smul g a
    have h2 : eB.toAddEquiv.symm (α g • e.toAddEquiv (eA.toAddEquiv a))
        = α g • eB.toAddEquiv.symm (e.toAddEquiv (eA.toAddEquiv a)) :=
      eB.symm.map_smul (α g) (e.toAddEquiv (eA.toAddEquiv a))
    simp only [AddEquiv.trans_apply, h1, e.map_smul]
    exact h2⟩

/-- 結論を存在文の形に落とす(`RecoverableAsAddModule` の形)。 -/
theorem exists_addEquiv {α : G ≃* G'} (e : SemilinearAddEquiv α A B) :
    ∃ φ : A ≃+ B, ∀ (g : G) (a : A), φ (g • a) = α g • φ a :=
  ⟨e.toAddEquiv, e.map_smul⟩

end SemilinearAddEquiv

/-! ## 2. ★★★抽象核 II —— 位相群上の局所定数関数(語彙ゼロ)

`LocallyConstant G M` は「`G` の位相群としての構造」と「加法群 `M`」だけから作られる。
したがって `G ≃ₜ* G'` はこれを**そのまま**運ぶ。これが本命題の「群論的」の内実である。 -/

section LocallyConstantSection

variable {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
variable {G' : Type*} [Group G'] [TopologicalSpace G'] [IsTopologicalGroup G']
variable {M M' : Type*} [AddCommGroup M] [AddCommGroup M']

/-- 左移動 `x ↦ g⁻¹ * x` を連続写像として。 -/
def leftTranslation (g : G) : C(G, G) :=
  ⟨fun x => g⁻¹ * x, continuous_const.mul continuous_id⟩

@[simp] theorem leftTranslation_apply (g x : G) : leftTranslation g x = g⁻¹ * x := rfl

/-- 平行移動 `(translate g f) x = f (g⁻¹ x)`。 -/
def translate (g : G) (f : LocallyConstant G M) : LocallyConstant G M :=
  LocallyConstant.comap (leftTranslation g) f

omit [AddCommGroup M] in
@[simp] theorem translate_apply (g : G) (f : LocallyConstant G M) (x : G) :
    translate g f x = f (g⁻¹ * x) := rfl

/-- 位相群 `G` 上の局所定数関数への `G` の**平行移動作用**。

★`scoped` にしてある:mathlib は `[SMul R Y]` から
`SMul R (LocallyConstant X Y)`(**値の側**への作用)を入れるので、
`R = G`・`X = G` のとき二つの作用が衝突しうる。 -/
scoped instance locallyConstantTranslation : DistribMulAction G (LocallyConstant G M) where
  smul := translate
  one_smul f := by ext x; show f (1⁻¹ * x) = f x; simp
  mul_smul g h f := by
    ext x
    show f ((g * h)⁻¹ * x) = f (h⁻¹ * (g⁻¹ * x))
    rw [mul_inv_rev, mul_assoc]
  smul_zero g := by
    ext x
    show (0 : LocallyConstant G M) (g⁻¹ * x) = (0 : LocallyConstant G M) x
    simp
  smul_add g f₁ f₂ := by
    ext x
    show (f₁ + f₂) (g⁻¹ * x) = f₁ (g⁻¹ * x) + f₂ (g⁻¹ * x)
    simp

@[simp] theorem locallyConstantTranslation_smul_apply (g : G) (f : LocallyConstant G M) (x : G) :
    (g • f) x = f (g⁻¹ * x) := rfl

/-- 係数の取り替え(加法同型の後合成)。 -/
def locallyConstantMapAddEquiv (e : M ≃+ M') :
    LocallyConstant G M ≃+ LocallyConstant G M' where
  toFun f := LocallyConstant.map e f
  invFun f := LocallyConstant.map e.symm f
  left_inv := by intro f; ext x; simp
  right_inv := by intro f; ext x; simp
  map_add' := by intro f₁ f₂; ext x; simp

omit [Group G] [IsTopologicalGroup G] in
@[simp] theorem locallyConstantMapAddEquiv_apply (e : M ≃+ M') (f : LocallyConstant G M) (x : G) :
    locallyConstantMapAddEquiv (G := G) e f x = e (f x) := rfl

/-- ★**係数の取り替えは平行移動作用と可換**。 -/
def locallyConstantCoeffCongr (e : M ≃+ M') :
    SemilinearAddEquiv (MulEquiv.refl G) (LocallyConstant G M) (LocallyConstant G M') :=
  ⟨locallyConstantMapAddEquiv e, by intro g f; ext x; rfl⟩

/-- `α : G ≃ₜ* G'` に沿った局所定数関数の引き戻し(加法同型)。 -/
def locallyConstantCongrLeft (α : G ≃ₜ* G') :
    LocallyConstant G M ≃+ LocallyConstant G' M where
  toFun f := LocallyConstant.comap
    (⟨α.toHomeomorph.symm, α.toHomeomorph.symm.continuous⟩ : C(G', G)) f
  invFun f := LocallyConstant.comap
    (⟨α.toHomeomorph, α.toHomeomorph.continuous⟩ : C(G, G')) f
  left_inv := by
    intro f; ext x
    show f (α.toHomeomorph.symm (α.toHomeomorph x)) = f x
    simp
  right_inv := by
    intro f; ext y
    show f (α.toHomeomorph (α.toHomeomorph.symm y)) = f y
    simp
  map_add' := by intro f₁ f₂; ext y; rfl

omit [IsTopologicalGroup G] [IsTopologicalGroup G'] in
@[simp] theorem locallyConstantCongrLeft_apply (α : G ≃ₜ* G') (f : LocallyConstant G M) (y : G') :
    locallyConstantCongrLeft (M := M) α f y = f (α.symm y) := rfl

/-- ★★★**本命題の抽象核**。位相群の同型 `α : G ≃ₜ* G'` は、局所定数関数のなす
加群の間の `α`-同変な加法同型を与える。

★分岐・付値・Galois の語彙が 1 語も出てこない。 -/
def locallyConstantSourceCongr (α : G ≃ₜ* G') :
    SemilinearAddEquiv α.toMulEquiv (LocallyConstant G M) (LocallyConstant G' M) :=
  ⟨locallyConstantCongrLeft α, by
    intro g f
    ext y
    show f (g⁻¹ * α.symm y) = f (α.symm ((α g)⁻¹ * y))
    have h : α.symm ((α g)⁻¹ * y) = g⁻¹ * α.symm y := by
      rw [map_mul, map_inv]
      simp
    rw [h]⟩

/-- 不変元はちょうど定数関数(★模型の係数が一意に決まることの根拠)。 -/
theorem smul_eq_self_iff_const (f : LocallyConstant G M) :
    (∀ g : G, g • f = f) ↔ ∃ m : M, f = LocallyConstant.const G m := by
  constructor
  · intro h
    refine ⟨f 1, ?_⟩
    ext x
    have hx := congrArg (fun u : LocallyConstant G M => u 1) (h x⁻¹)
    simp only [locallyConstantTranslation_smul_apply, inv_inv, mul_one] at hx
    show f x = f 1
    exact hx
  · rintro ⟨m, rfl⟩
    intro g
    ext x
    rfl

end LocallyConstantSection

/-! ## 3. 具体層 —— `K̄` の群論的模型 -/

variable {p : ℕ} [Fact p.Prime]

/-- `Γ_K` の位相群としての構造と自然数 `d` だけから作れる加群 —— `ℚ_p^d` に値をとる
`Γ_K` 上の局所定数関数。★体 `K` は台にも値にも現れない。 -/
abbrev smoothModelType (K : PAdicLocalField p) (d : ℕ) : Type :=
  LocallyConstant K.absGal (Fin d → ℚ_[p])

/-- **群論的模型を持つ**:`K̄` が `Γ_K`-同変に `smoothModelType K d` と同型。

★これは「無限次 Galois 拡大に対する正規底定理」に `K ≅ ℚ_p^d` を合わせたもの。 -/
def SmoothModel (K : PAdicLocalField p) (d : ℕ) : Prop :=
  Nonempty (SemilinearAddEquiv (MulEquiv.refl K.absGal) K.closure (smoothModelType K d))

/-- 原典の言葉に忠実な形 —— `K̄ ≅ C^∞(Γ_K, K)`(無限次版の正規底定理そのもの)。 -/
def SmoothModelCarrier (K : PAdicLocalField p) : Prop :=
  Nonempty (SemilinearAddEquiv (MulEquiv.refl K.absGal) K.closure
    (LocallyConstant K.absGal K.carrier))

/-- 忠実な形から `ℚ_p^d` 版へ —— `K ≅ ℚ_p^{[K:ℚ_p]}`(有限次元線型代数)を係数に入れるだけ。 -/
theorem smoothModel_of_carrier (K : PAdicLocalField p) (h : SmoothModelCarrier K) :
    SmoothModel K (Module.finrank ℚ_[p] K.carrier) := by
  obtain ⟨e⟩ := h
  refine ⟨SemilinearAddEquiv.ofEq ?_
    (e.trans (locallyConstantCoeffCongr
      (G := K.absGal) (Module.finBasis ℚ_[p] K.carrier).equivFun.toAddEquiv))⟩
  ext g
  rfl

/-- ★★★**[pGC] Proposition 2.1 への還元**。

各 `K` について `K̄` が群論的模型 `LocallyConstant Γ_K (ℚ_p^{[K:ℚ_p]})` と
`Γ_K`-同変に同型でありさえすれば、`K̄` は `Γ_K` から群論的に回復される。

`[K:ℚ_p] = [K':ℚ_p]` は `finrank_eq_of_absGal_equiv`(`Found/PGC/DegreeTransport.lean`、
`sorry` 無し)から出る —— ★ここで §1 の結果が効いている。 -/
theorem prop_2_1_of_smoothModel
    (h : ∀ K : PAdicLocalField p, SmoothModel K (Module.finrank ℚ_[p] K.carrier)) :
    RecoverableAsAddModule (p := p) (fun K => K.closure) := by
  intro K K' α
  have hd : Module.finrank ℚ_[p] K.carrier = Module.finrank ℚ_[p] K'.carrier :=
    finrank_eq_of_absGal_equiv K K' α
  obtain ⟨eK⟩ := hd ▸ h K
  obtain ⟨eK'⟩ := h K'
  exact (SemilinearAddEquiv.sandwich eK (locallyConstantSourceCongr α) eK').exists_addEquiv

/-- 忠実な形からの版。 -/
theorem prop_2_1_of_smoothModelCarrier
    (h : ∀ K : PAdicLocalField p, SmoothModelCarrier K) :
    RecoverableAsAddModule (p := p) (fun K => K.closure) :=
  prop_2_1_of_smoothModel (fun K => smoothModel_of_carrier K (h K))

/-! ## 4. ★非空虚性 —— 還元先は空虚ではない -/

/-- ★**局所定数関数の族そのものは、無条件に群論的に回復される**(`sorry` 無し)。

`Obj K := LocallyConstant Γ_K M` に対して `RecoverableAsAddModule` が成り立つ。
これは `SmoothModel` という形の仮説が空虚でないこと(そういう `Γ_K`-加群が
実際に存在すること)の witness である。 -/
theorem recoverableAsAddModule_locallyConstant (M : Type) [AddCommGroup M] :
    RecoverableAsAddModule (p := p) (fun K => LocallyConstant K.absGal M) := by
  intro K K' α
  exact (locallyConstantSourceCongr (M := M) α).exists_addEquiv

/-- 同上を模型の型で。 -/
theorem recoverableAsAddModule_smoothModelType (d : ℕ) :
    RecoverableAsAddModule (p := p) (fun K => smoothModelType K d) :=
  recoverableAsAddModule_locallyConstant (p := p) (Fin d → ℚ_[p])

/-! ## 5. ★★有限次 Galois 拡大での模型 —— mathlib の正規底定理の**同変形**

mathlib の `IsGalois.normalBasis` は「Galois 群の軌道になっている基底が存在する」
までしか言わない。本節はそれを `Gal(E/F)`-**加群の同型**に上げる:

  `E ≅ LocallyConstant Gal(E/F) F`  (`Gal(E/F)`-同変)

有限次拡大では `Gal(E/F)` は離散(`krullTopology_discreteTopology_of_finiteDimensional`)
なので `LocallyConstant Gal(E/F) F` は「すべての関数」である。
★これは `SmoothModelCarrier` の**有限次の段**そのものである。 -/

section FiniteLevel

variable (F E : Type*) [Field F] [Field E] [Algebra F E] [FiniteDimensional F E] [IsGalois F E]

/-- 離散位相の空間上では任意の関数が局所定数。 -/
def ofDiscreteFun {X Y : Type*} [TopologicalSpace X] [DiscreteTopology X] (f : X → Y) :
    LocallyConstant X Y := ⟨f, fun _ => isOpen_discrete _⟩

@[simp] theorem ofDiscreteFun_apply {X Y : Type*} [TopologicalSpace X] [DiscreteTopology X]
    (f : X → Y) (x : X) : ofDiscreteFun f x = f x := rfl

/-- 正規底は Galois 群の作用で**置換される**:`σ (b τ) = b (σ τ)`。

`IsGalois.normalBasis_apply`(`b e = e (b 1)`)の直接の帰結。 -/
theorem normalBasis_apply_smul (σ τ : E ≃ₐ[F] E) :
    σ (IsGalois.normalBasis F E τ) = IsGalois.normalBasis F E (σ * τ) := by
  rw [IsGalois.normalBasis_apply τ, IsGalois.normalBasis_apply (σ * τ), AlgEquiv.mul_apply]

/-- 正規底での座標表示は平行移動と同変:`[σ x]_ρ = [x]_{σ⁻¹ρ}`。 -/
theorem normalBasis_equivFun_smul (σ : E ≃ₐ[F] E) (x : E) (ρ : E ≃ₐ[F] E) :
    (IsGalois.normalBasis F E).equivFun (σ x) ρ
      = (IsGalois.normalBasis F E).equivFun x (σ⁻¹ * ρ) := by
  classical
  set b := IsGalois.normalBasis F E with hb
  have h : (b.equivFun.toLinearMap.comp σ.toLinearEquiv.toLinearMap)
      = (LinearMap.funLeft F F (fun r : E ≃ₐ[F] E => σ⁻¹ * r)).comp b.equivFun.toLinearMap := by
    refine b.ext fun τ => ?_
    funext ρ'
    simp only [LinearMap.comp_apply, LinearEquiv.coe_coe, AlgEquiv.toLinearEquiv_apply,
      LinearMap.funLeft_apply]
    rw [hb, normalBasis_apply_smul, ← hb, b.equivFun_self, b.equivFun_self]
    exact if_congr eq_inv_mul_iff_mul_eq.symm rfl rfl
  exact congrFun (congrArg (fun m : E →ₗ[F] ((E ≃ₐ[F] E) → F) => m x) h) ρ

/-- ★★★**正規底定理の同変形**。有限次 Galois 拡大 `E/F` について、
`E` は `Gal(E/F)`-加群として `Gal(E/F)` 上の `F`-値(局所定数)関数と同型。

mathlib が与えるのは基底 `IsGalois.normalBasis` までで、
`Gal(E/F)`-同変性はここで初めて付く。 -/
noncomputable def finiteNormalBasisModel :
    SemilinearAddEquiv (MulEquiv.refl (E ≃ₐ[F] E)) E (LocallyConstant (E ≃ₐ[F] E) F) :=
  ⟨{ toFun := fun x => ofDiscreteFun ((IsGalois.normalBasis F E).equivFun x)
     invFun := fun f => (IsGalois.normalBasis F E).equivFun.symm (f : (E ≃ₐ[F] E) → F)
     left_inv := by
       intro x
       show (IsGalois.normalBasis F E).equivFun.symm ((IsGalois.normalBasis F E).equivFun x) = x
       simp
     right_inv := by
       intro f
       ext ρ
       show (IsGalois.normalBasis F E).equivFun
         ((IsGalois.normalBasis F E).equivFun.symm (f : (E ≃ₐ[F] E) → F)) ρ = f ρ
       rw [LinearEquiv.apply_symm_apply]
     map_add' := by
       intro x y
       ext ρ
       show (IsGalois.normalBasis F E).equivFun (x + y) ρ
         = (IsGalois.normalBasis F E).equivFun x ρ + (IsGalois.normalBasis F E).equivFun y ρ
       simp },
    by
      intro σ x
      ext ρ
      show (IsGalois.normalBasis F E).equivFun (σ x) ρ
        = (IsGalois.normalBasis F E).equivFun x (σ⁻¹ * ρ)
      exact normalBasis_equivFun_smul F E σ x ρ⟩

end FiniteLevel

end ABC3.Found.PGC

