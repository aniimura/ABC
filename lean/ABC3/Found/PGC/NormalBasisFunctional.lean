import ABC3.Found.PGC.SmoothModelTransport
import Mathlib.FieldTheory.Galois.Profinite

/-!
# 汎関数による正規底模型 —— `K̄ ≅ C^∞(Γ_K, K)` への**塔をまたぐ**道具立て

本ファイルは `Found/PGC/SmoothModelTransport.lean` が残した仮説

  `SmoothModelCarrier K` :`K̄ ≅ LocallyConstant Γ_K K`(`Γ_K`-同変)

を、**有限次の水準に完全に分解する**。要は次の 1 本の写像である:

  `Φ_λ : A → (G → M)`,  `Φ_λ (a) (g) = λ (g⁻¹ • a)`   (`λ : A →+ M` は加法的汎関数)

## ★★★なぜ「正規底の生成元」ではなく「汎関数」なのか(前段の設計変更)

前段(`SmoothModelTransport.lean` の docstring)は、無限次への持ち上げを

  「有限次 Galois `L ⊆ L'` ごとに正規底生成元 `θ_L` を選び、
   `θ_L = Tr_{L'/L} θ_{L'}` という**両立系**を作る」

と描いていた。この道では遷移の可換性(下の落とし穴)を毎回確かめる必要があり、
さらに `K[G'] ↠ K[G]` の**単元の持ち上げ**(Maschke)を要する。

★**汎関数の側に移すと、遷移の可換性が 1 行になる。**
`ι : B →+ A` が `π : G →* Q` に沿って同変なら、`μ := λ ∘ ι` について

  `Φ_λ (ι b) = Φ_μ (b) ∘ π`      (`OrbitModel.orbitFun_comp`、証明 3 行)

が**無条件に**成り立つ。すなわち

* 体の側の遷移は**包含** `L ↪ L'`、
* 模型の側の遷移は `Gal(L'/K) ↠ Gal(L/K)` に沿った**引き戻し**、

という向きの食い違い(本持ち場で最初に検算すべきと指摘された点)は、
「1 つの `λ` を制限する」という形にすると**自動的に整合する** ——
`λ` の制限は関手的だから、両立系を選ぶ必要がそもそも無い。

★さらに**降下**(`OrbitModel.bijective_descent`)も Maschke も跡写像も使わずに出る:
`Φ_λ` が全体で全単射なら、各水準 `A^{ker π} = range ι` の上で `Φ_{λ∘ι}` は全単射である。
(前段が「`K[G']^H ≅ K[G]` の下で基底元は基底元に落ちる」と述べていた事実の、
群環を経由しない証明。)

## この分解で**残る**もの(★正直な記録)

本ファイルの結論は**同値**である(`smoothModelCarrier_iff_hasCoherentFunctional`):

  `SmoothModelCarrier K`  ⟺  `HasCoherentFunctional K`
  (`HasCoherentFunctional K` := ∃ λ : K̄ →+ K で、**すべての**有限次 Galois `L/K` に対し
   `Φ_{λ|_L} : L → Maps(Gal(L/K), K)` が全単射)

★**逆向きも証明した**ので、この仮説は強すぎも弱すぎもしない ——
残っている数学は**ちょうどこれだけ**である。
(逆向きの鍵は Frobenius 相互律:`Γ_K`-同変な加法同型 `K̄ ≅ C^∞(Γ_K,K)` は
**必ず** `Φ_λ` の形をしている、`λ := (1 での評価) ∘ E`。)

各水準**単独**では条件は必ず満たせる(`exists_levelwise_functional`、mathlib の
正規底定理)。難しいのは「**1 つの** `λ` で全水準を同時に」であり、そこには

1. **Maschke による上げ**:`λ|_L` が良いとき、`L ⊆ L'` に良いまま延長できること。
   `K[G'] = K[G] × A`(char 0 で半単純)の中心冪等元
   `e = (1/|N|) Σ_{n ∈ N} n` を使う。★その**環論だけの核**は本ファイルに在る
   (`OrbitModel.exists_unit_of_central_idempotent`、群も体も出てこない)。
2. **可算性**:`K̄/K` の有限次部分拡大が可算個しかないこと
   (各次数について有限個 —— Krasner の補題 + 係数の compactness)。
   ★有向集合上の全射逆系は極限が空になりうるので、可算共終列がどうしても要る。
   ★2026-09-08 実測:**Krasner の補題そのものは mathlib に在る**
   (`grep -n -i "krasner" .cache/mathlib-index.txt` →
   `IsKrasner` / `IsKrasner.krasner` / `IsKrasner.of_completeSpace`、
   `Mathlib/Analysis/Normed/Field/Krasner.lean`)。
   ★無いのは**その系**「与えられた次数の拡大は有限個」「有限次部分拡大は可算」で、
   `grep -n "finite_extensions\|countable.*IntermediateField" .cache/mathlib-index.txt`
   は 0 件だった。本木にも無い。★これが新しい節点である。

の 2 つが要る。本ファイルはそこまでは埋めない。

★★**2026-09-08 追記(後続の波):この 2 つは両方とも埋まった。**

* (1) Maschke 上げ:`Found/PGC/CoherentFunctional.lean::OrbitModel.exists_bijective_extension`。
  ★**群環との同一視は要らなかった** —— 同変な収縮 `ρ` があれば
  `Ψ := ι ∘ φ ∘ ρ + (1 − ι ∘ ρ)` が同変自己同型になる、という初等的な観察で足りる。
  ★上の `exists_unit_of_central_idempotent` は**使っていない**。
* (2) 可算性:`Found/PGC/CountableGenerators.lean::hasCountableGenerators`。
  ★**ℚ 上代数的な元は要らなかった** —— 係数を `K` の可算稠密部分集合から取れば
  Krasner の補題(`IsKrasner.krasner`、mathlib)がそのまま効く。

したがって `HasCoherentFunctional K` は無条件に成り立ち
(`Found/PGC/CountableGenerators.lean::hasCoherentFunctional`)、
[pGC] Proposition 2.1 は閉じた(同ファイル `prop_2_1`)。

## ★真偽の検算(本持ち場で最初に指示されたこと)

主張は**真**である。検算した点:

* **遷移の向き**:体の側は包含 `L ↪ L'`、模型の側は `Gal(L'/K) ↠ Gal(L/K)` に沿った
  引き戻し —— この 2 つが可換になるかが要と指摘された。★汎関数で書くと
  `orbitFun_comp` が示すとおり**無条件に可換**である(`λ` を制限するだけだから)。
  正規底生成元で書いたときの整合条件は `θ_L = Tr_{L'/L} θ_{L'}` で、
  実際 `Tr θ'` は `δ_1 ∈ Maps(G',K)` の `H` 上の和 = `δ_1` の引き戻しに対応する。
  ★どちらの記述も同じものだが、汎関数の側では**選択が要らない**。
* **`K[G]`(左移動)と `Maps(G,K)` の同一視**:`Σ a_g g ↦ (x ↦ a_x)` が同変
  (`g · Σ a_h h = Σ a_{g⁻¹x} x`、`(g·f)(x) = f(g⁻¹x)`)。逆元の位置は
  `orbitFun`(`Φ_λ(a)(g) = λ(g⁻¹ • a)`)に埋め込んであり、
  `orbitFun_smul` が同変性を保証する。
* **Maschke で単元が持ち上がるか**:「分裂する」より強い主張だと指摘されたとおりで、
  必要なのは中心冪等元 `e` による**角(corner)分解**である。
  `exists_unit_of_central_idempotent` はまさにその形で、
  持ち上げは `u ↦ a·e + (1 - e)` (`a` は任意の持ち上げ)。

## 内容

§1 抽象核 I:`Φ_λ` そのもの(群・加法群だけ。分岐・付値・Galois の語彙ゼロ)。
§2 抽象核 II:遷移の可換性と降下(同上)。局所定数版と単元持ち上げ核も同節。
§3 抽象核 III:コンパクト位相群上の局所定数関数の一様性(同上)。
§4 具体層(有限次):mathlib の `IsGalois.normalBasis` から
   `Φ_{coord 1}` が全単射であること —— ★水準ごとの条件の**非空虚性**。
§5 具体層(無限次):`Γ_K` の水準への分解、`SmoothModelCarrier` との**同値**、
   および [pGC] Proposition 2.1 への配線(`prop_2_1_of_hasCoherentFunctional`)。

★**在庫の測り方の記録**(2026-09-08):`IsGalois.normalBasis` は
`.cache/mathlib-index.txt` に**載っていない**(同ファイル `FieldTheory/Galois/NormalBasis.lean`
からは 2 つの `exists_linearIndependent_...` しか索引されていない)。
索引に無いことは不在の証拠にならない —— `#check @IsGalois.normalBasis` で在ることを確かめた。
-/

namespace ABC3.Found.PGC

/-! ## 1. ★★★抽象核 I —— 汎関数の軌道写像(語彙ゼロ) -/

namespace OrbitModel

variable {G : Type*} [Group G] {A : Type*} [AddCommGroup A] [DistribMulAction G A]
variable {M : Type*} [AddCommGroup M]

/-- **汎関数 `λ` の軌道写像** `Φ_λ(a) : g ↦ λ (g⁻¹ • a)`。

★これは `Hom_G(A, Maps(G,M)) ≅ Hom(A,M)`(Frobenius 相互律)の右辺から左辺への写像であり、
`G`-同変性は**作り付け**である(`orbitFun_smul`)。 -/
def orbitFun (lam : A →+ M) (a : A) : G → M := fun g => lam (g⁻¹ • a)

@[simp] theorem orbitFun_apply (lam : A →+ M) (a : A) (g : G) :
    orbitFun lam a g = lam (g⁻¹ • a) := rfl

/-- ★同変性(平行移動作用 `(h • f) g = f (h⁻¹ * g)` に対して)。 -/
theorem orbitFun_smul (lam : A →+ M) (h : G) (a : A) (g : G) :
    orbitFun lam (h • a) g = orbitFun lam a (h⁻¹ * g) := by
  simp [orbitFun, smul_smul, mul_inv_rev]

/-- 加法群準同型としての `Φ_λ`。 -/
def orbitHom (lam : A →+ M) : A →+ (G → M) where
  toFun := orbitFun lam
  map_zero' := by funext g; simp [orbitFun]
  map_add' a b := by funext g; simp [orbitFun]

@[simp] theorem orbitHom_apply (lam : A →+ M) (a : A) (g : G) :
    orbitHom (G := G) lam a g = lam (g⁻¹ • a) := rfl

/-- 単射性は「`λ` の `G`-軌道が `A` を分離する」ことと同値。 -/
theorem orbitHom_injective_iff (lam : A →+ M) :
    Function.Injective (orbitHom (G := G) lam) ↔ ∀ a : A, (∀ g : G, lam (g • a) = 0) → a = 0 := by
  rw [injective_iff_map_eq_zero]
  constructor
  · intro h a ha
    refine h a ?_
    ext g
    simpa using ha g⁻¹
  · intro h a ha
    refine h a fun g => ?_
    have := congrFun (congrArg (fun u : G → M => u) ha) g⁻¹
    simpa using this

section Topological

variable [TopologicalSpace G] [IsTopologicalGroup G]

/-- 固定部分群が開なら軌道写像は局所定数(右剰余類 `S·g` の上で定数)。 -/
theorem isLocallyConstant_orbitFun (lam : A →+ M) (a : A)
    (h : IsOpen ((MulAction.stabilizer G a : Subgroup G) : Set G)) :
    IsLocallyConstant (orbitFun (G := G) lam a) := by
  rw [IsLocallyConstant.iff_exists_open]
  intro g
  refine ⟨(fun x => x * g) '' ((MulAction.stabilizer G a : Subgroup G) : Set G), ?_, ?_, ?_⟩
  · exact (Homeomorph.mulRight g).isOpenMap _ h
  · exact ⟨1, (MulAction.stabilizer G a).one_mem, by simp⟩
  · rintro _ ⟨s, hs, rfl⟩
    have hsa : s⁻¹ • a = a := (MulAction.stabilizer G a).inv_mem hs
    show lam ((s * g)⁻¹ • a) = lam (g⁻¹ • a)
    rw [mul_inv_rev, ← smul_smul, hsa]

/-- 局所定数関数の加群への `Φ_λ`。 -/
def orbitLC (lam : A →+ M) (hlc : ∀ a : A, IsLocallyConstant (orbitFun (G := G) lam a)) :
    A →+ LocallyConstant G M where
  toFun a := ⟨orbitFun lam a, hlc a⟩
  map_zero' := by ext g; simp
  map_add' a b := by ext g; simp

omit [IsTopologicalGroup G] in
@[simp] theorem orbitLC_apply (lam : A →+ M)
    (hlc : ∀ a : A, IsLocallyConstant (orbitFun (G := G) lam a)) (a : A) (g : G) :
    orbitLC lam hlc a g = lam (g⁻¹ • a) := rfl

/-- ★`Φ_λ` は作り付けで `G`-同変。 -/
theorem orbitLC_map_smul (lam : A →+ M)
    (hlc : ∀ a : A, IsLocallyConstant (orbitFun (G := G) lam a)) (h : G) (a : A) :
    orbitLC lam hlc (h • a) = h • orbitLC lam hlc a := by
  ext g
  show orbitFun lam (h • a) g = orbitFun lam a (h⁻¹ * g)
  exact orbitFun_smul lam h a g

/-- ★★`Φ_λ` が全単射なら、それが求める `SemilinearAddEquiv` である。 -/
noncomputable def semilinearOfBijective (lam : A →+ M)
    (hlc : ∀ a : A, IsLocallyConstant (orbitFun (G := G) lam a))
    (hbij : Function.Bijective (orbitLC lam hlc)) :
    SemilinearAddEquiv (MulEquiv.refl G) A (LocallyConstant G M) :=
  ⟨AddEquiv.ofBijective (orbitLC lam hlc) hbij, fun g a => orbitLC_map_smul lam hlc g a⟩

end Topological

/-! ## 2. ★★★抽象核 II —— 遷移の可換性と降下(語彙ゼロ) -/

section Levels

variable {Q : Type*} [Group Q] {B : Type*} [AddCommGroup B] [DistribMulAction Q B]

/-- ★★★**遷移の可換性**。`ι : B →+ A` が全射 `π : G →* Q` に沿って同変なら

  `Φ_λ ∘ ι = (π に沿った引き戻し) ∘ Φ_{λ∘ι}`。

★本持ち場で最初に検算すべきと言われた点そのもの:体の側の遷移(包含 `ι`)と
模型の側の遷移(`π` に沿った引き戻し)が、汎関数の制限の下で**自動的に**可換になる。 -/
theorem orbitFun_comp (pi : G →* Q) (iota : B →+ A)
    (hiota : ∀ (g : G) (b : B), iota (pi g • b) = g • iota b)
    (lam : A →+ M) (b : B) :
    orbitFun (G := G) lam (iota b) = (orbitFun (G := Q) (lam.comp iota) b) ∘ pi := by
  funext g
  show lam (g⁻¹ • iota b) = lam (iota ((pi g)⁻¹ • b))
  rw [← map_inv pi g, hiota]

/-- 水準を経由する関数は `range ι` の元で実現される(全射性の材料)。 -/
theorem exists_orbitFun_eq_of_level (pi : G →* Q) (iota : B →+ A)
    (hiota : ∀ (g : G) (b : B), iota (pi g • b) = g • iota b) (lam : A →+ M)
    (hlev : Function.Surjective (orbitHom (G := Q) (lam.comp iota))) (f : Q → M) :
    ∃ b : B, orbitFun (G := G) lam (iota b) = f ∘ pi := by
  obtain ⟨b, hb⟩ := hlev f
  refine ⟨b, ?_⟩
  rw [orbitFun_comp pi iota hiota lam b]
  exact congrArg (fun u : Q → M => u ∘ pi) hb

/-- 水準での単射性は `range ι` の元を分離する(単射性の材料)。 -/
theorem eq_of_orbitFun_eq_of_level (pi : G →* Q) (hpi : Function.Surjective pi) (iota : B →+ A)
    (hiota : ∀ (g : G) (b : B), iota (pi g • b) = g • iota b) (lam : A →+ M)
    (hlev : Function.Injective (orbitHom (G := Q) (lam.comp iota))) (b₁ b₂ : B)
    (h : orbitFun (G := G) lam (iota b₁) = orbitFun (G := G) lam (iota b₂)) : b₁ = b₂ := by
  refine hlev ?_
  funext q
  obtain ⟨g, rfl⟩ := hpi q
  have h1 := congrFun h g
  rw [orbitFun_comp pi iota hiota lam b₁, orbitFun_comp pi iota hiota lam b₂] at h1
  exact h1

/-- ★★★**降下**。全体で全単射なら各水準でも全単射。

★Maschke も跡写像も群環も使わない —— 使うのは
`A^{ker π} ⊆ range ι`(Galois 降下)と `π` の全射性だけである。
前段が「跡写像 `Tr_{L'/L}` で基底元が基底元に落ちる」と述べていた事実の代替。 -/
theorem bijective_descent (pi : G →* Q) (hpi : Function.Surjective pi) (iota : B →+ A)
    (hiotainj : Function.Injective iota)
    (hiota : ∀ (g : G) (b : B), iota (pi g • b) = g • iota b)
    (hfix : ∀ a : A, (∀ n : G, pi n = 1 → n • a = a) → a ∈ Set.range iota)
    (lam : A →+ M) (hglob : Function.Bijective (orbitHom (G := G) lam)) :
    Function.Bijective (orbitHom (G := Q) (lam.comp iota)) := by
  constructor
  · intro b₁ b₂ h
    refine hiotainj (hglob.1 ?_)
    funext g
    show orbitFun (G := G) lam (iota b₁) g = orbitFun (G := G) lam (iota b₂) g
    rw [orbitFun_comp pi iota hiota lam b₁, orbitFun_comp pi iota hiota lam b₂]
    exact congrFun (congrArg (fun u : Q → M => u ∘ pi) h) g
  · intro f
    obtain ⟨a, ha⟩ := hglob.2 (f ∘ pi)
    have e1 : orbitFun (G := G) lam a = f ∘ pi := ha
    have hfixed : ∀ n : G, pi n = 1 → n • a = a := by
      intro n hn
      refine hglob.1 ?_
      funext g
      show orbitFun (G := G) lam (n • a) g = orbitFun (G := G) lam a g
      rw [orbitFun_smul lam n a g, e1]
      show f (pi (n⁻¹ * g)) = f (pi g)
      rw [map_mul, map_inv, hn]
      simp
    obtain ⟨b, rfl⟩ := hfix a hfixed
    refine ⟨b, ?_⟩
    funext q
    obtain ⟨g, rfl⟩ := hpi q
    show orbitFun (G := Q) (lam.comp iota) b (pi g) = f (pi g)
    have h2 : orbitFun (G := Q) (lam.comp iota) b (pi g) = orbitFun (G := G) lam (iota b) g :=
      congrFun (orbitFun_comp pi iota hiota lam b).symm g
    rw [h2, e1]
    rfl

end Levels

/-! ### 2.1 局所定数版の降下(★大域が `LocallyConstant` に限られている場合) -/

section LevelsTopological

variable [TopologicalSpace G] [IsTopologicalGroup G]
variable {Q : Type*} [Group Q] {B : Type*} [AddCommGroup B] [DistribMulAction Q B]

omit [AddCommGroup M] in
/-- 核が開なら、水準を経由する関数はすべて局所定数。 -/
theorem isLocallyConstant_comp_of_isOpen_ker (pi : G →* Q)
    (h : IsOpen ((pi.ker : Subgroup G) : Set G)) (f : Q → M) :
    IsLocallyConstant (f ∘ pi) := by
  rw [IsLocallyConstant.iff_exists_open]
  intro g
  refine ⟨(fun x => g * x) '' ((pi.ker : Subgroup G) : Set G),
    (Homeomorph.mulLeft g).isOpenMap _ h, ⟨1, one_mem _, by simp⟩, ?_⟩
  rintro _ ⟨s, hs, rfl⟩
  show f (pi (g * s)) = f (pi g)
  rw [map_mul, MonoidHom.mem_ker.1 hs, mul_one]

omit [IsTopologicalGroup G] in
/-- ★★`bijective_descent` の**局所定数版**。大域の全単射性が
`orbitLC`(局所定数関数への写像)についてのものであっても、水準では
`orbitHom`(すべての関数への写像)が全単射になる。 -/
theorem bijective_descent_lc (pi : G →* Q) (hpi : Function.Surjective pi) (iota : B →+ A)
    (hiotainj : Function.Injective iota)
    (hiota : ∀ (g : G) (b : B), iota (pi g • b) = g • iota b)
    (hfix : ∀ a : A, (∀ n : G, pi n = 1 → n • a = a) → a ∈ Set.range iota)
    (hlcpi : ∀ f : Q → M, IsLocallyConstant (f ∘ pi))
    (lam : A →+ M) (hlc : ∀ a : A, IsLocallyConstant (orbitFun (G := G) lam a))
    (hglob : Function.Bijective (orbitLC lam hlc)) :
    Function.Bijective (orbitHom (G := Q) (lam.comp iota)) := by
  have hinj : Function.Injective (orbitFun (G := G) lam) := by
    intro x y hxy
    exact hglob.1 (by ext g; exact congrFun hxy g)
  constructor
  · intro b₁ b₂ h
    refine hiotainj (hinj ?_)
    funext g
    rw [orbitFun_comp pi iota hiota lam b₁, orbitFun_comp pi iota hiota lam b₂]
    exact congrFun (congrArg (fun u : Q → M => u ∘ pi) h) g
  · intro f
    obtain ⟨a, ha⟩ := hglob.2 ⟨f ∘ pi, hlcpi f⟩
    have e1 : orbitFun (G := G) lam a = f ∘ pi := by
      funext g
      exact congrFun (congrArg (fun u : LocallyConstant G M => (u : G → M)) ha) g
    have hfixed : ∀ n : G, pi n = 1 → n • a = a := by
      intro n hn
      refine hinj ?_
      funext g
      rw [orbitFun_smul lam n a g, e1]
      show f (pi (n⁻¹ * g)) = f (pi g)
      rw [map_mul, map_inv, hn]
      simp
    obtain ⟨b, rfl⟩ := hfix a hfixed
    refine ⟨b, ?_⟩
    funext q
    obtain ⟨g, rfl⟩ := hpi q
    show orbitFun (G := Q) (lam.comp iota) b (pi g) = f (pi g)
    have h2 : orbitFun (G := Q) (lam.comp iota) b (pi g) = orbitFun (G := G) lam (iota b) g :=
      congrFun (orbitFun_comp pi iota hiota lam b).symm g
    rw [h2, e1]
    rfl

end LevelsTopological

/-! ### 2.2 ★★★抽象核 —— 中心冪等元に沿った単元の持ち上げ(Maschke の受け皿)

`K[G'] ↠ K[G]` の単元が持ち上がる、という主張の**環論だけの核**。
`e := (1/|N|) Σ_{n ∈ N} n` が中心冪等元で `f e = 1`、`ker f` が右から `e` を消すとき、
`u ↦ (a e + (1 - e))` が持ち上げになる。★群も体も出てこない。 -/

section UnitLift

variable {R S : Type*} [Ring R] [Ring S]

theorem exists_unit_of_central_idempotent (f : R →+* S) (hf : Function.Surjective f) {e : R}
    (he : e * e = e) (hc : ∀ a : R, e * a = a * e) (hfe : f e = 1)
    (hker : ∀ x : R, f x = 0 → x * e = 0) (u : Sˣ) : ∃ v : Rˣ, f v = u := by
  obtain ⟨a, ha⟩ := hf u
  obtain ⟨b, hb⟩ := hf (u⁻¹ : Sˣ)
  have habe : ∀ x y : R, f (x * y) = 1 → x * y * e = e := by
    intro x y hxy
    have h0 : f (x * y - 1) = 0 := by rw [map_sub, hxy, map_one, sub_self]
    have h1 := hker _ h0
    rw [sub_mul, one_mul, sub_eq_zero] at h1
    exact h1
  have expand : ∀ x y : R, x * y * e = e → (x * e + (1 - e)) * (y * e + (1 - e)) = 1 := by
    intro x y hxy
    have h1 : x * e * (y * e) = e := by
      have hey : e * (y * e) = y * (e * e) := by rw [← mul_assoc, hc y, mul_assoc]
      rw [mul_assoc, hey, he, ← mul_assoc]
      exact hxy
    have h2 : x * e * (1 - e) = 0 := by
      rw [mul_sub, mul_one, mul_assoc, he, sub_self]
    have h3 : (1 - e) * (y * e) = 0 := by
      have hey : e * (y * e) = y * e := by rw [← mul_assoc, hc y, mul_assoc, he]
      rw [sub_mul, one_mul, hey, sub_self]
    have h4 : (1 - e) * (1 - e) = 1 - e := by
      have hx : (1 - e) * (1 - e) = 1 - e - e + e * e := by noncomm_ring
      rw [hx, he]
      abel
    have hexp : (x * e + (1 - e)) * (y * e + (1 - e))
        = x * e * (y * e) + x * e * (1 - e) + ((1 - e) * (y * e) + (1 - e) * (1 - e)) := by
      noncomm_ring
    rw [hexp, h1, h2, h3, h4]
    abel
  refine ⟨⟨a * e + (1 - e), b * e + (1 - e), expand a b (habe a b ?_), expand b a (habe b a ?_)⟩, ?_⟩
  · rw [map_mul, ha, hb]
    exact u.mul_inv
  · rw [map_mul, hb, ha]
    exact u.inv_mul
  · show f (a * e + (1 - e)) = (u : S)
    rw [map_add, map_mul, ha, hfe, map_sub, map_one, hfe, mul_one, sub_self, add_zero]

end UnitLift

/-! ## 3. ★★★抽象核 III —— コンパクト位相群上の局所定数関数(語彙ゼロ) -/

section Compact

variable {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G] [CompactSpace G]
variable {M : Type*}

/-- **コンパクト位相群上の局所定数関数は 1 のある近傍で右一様**:
`∃ U ∈ 𝓝 1, ∀ g u ∈ U, f (g * u) = f g`。

★これが「局所定数関数は有限商を経由する」の中身である(profinite の圏論的極限を使わない)。 -/
theorem exists_nhds_one_forall_mul (f : LocallyConstant G M) :
    ∃ U ∈ nhds (1 : G), ∀ (g u : G), u ∈ U → f (g * u) = f g := by
  classical
  have key : ∀ g : G, ∃ V : Set G, IsOpen V ∧ (1 : G) ∈ V ∧
      ∀ v ∈ V, ∀ w ∈ V, f (g * (v * w)) = f g := by
    intro g
    have hop : IsOpen {v : G | f (g * v) = f g} :=
      (f.isLocallyConstant.isOpen_fiber (f g)).preimage (continuous_const.mul continuous_id)
    have hmem : {v : G | f (g * v) = f g} ∈ nhds (1 : G) := hop.mem_nhds (by simp)
    obtain ⟨V, hVopen, hV1, hVsub⟩ := exists_open_nhds_one_mul_subset hmem
    exact ⟨V, hVopen, hV1, fun v hv w hw => hVsub (Set.mul_mem_mul hv hw)⟩
  choose V hVopen hV1 hVspec using key
  have hcov : (Set.univ : Set G) ⊆ ⋃ g : G, (fun x => g * x) '' (V g) := by
    intro g _
    exact Set.mem_iUnion.2 ⟨g, ⟨1, hV1 g, by simp⟩⟩
  obtain ⟨t, ht⟩ := isCompact_univ.elim_finite_subcover (fun g : G => (fun x => g * x) '' (V g))
    (fun g => (isOpenMap_mul_left g) _ (hVopen g)) hcov
  refine ⟨⋂ g ∈ t, V g, ?_, ?_⟩
  · exact (Filter.biInter_finset_mem t).2 fun i _ => (hVopen i).mem_nhds (hV1 i)
  · intro g u hu
    obtain ⟨i, hit, v, hv, hgv⟩ := Set.mem_iUnion₂.1 (ht (Set.mem_univ g))
    have hui : u ∈ V i := Set.mem_iInter₂.1 hu i hit
    have h1 : f (g * u) = f i := by
      rw [← hgv, mul_assoc]
      exact hVspec i v hv u hui
    have h2 : f g = f i := by
      rw [← hgv]
      have := hVspec i v hv 1 (hV1 i)
      rwa [mul_one] at this
    rw [h1, h2]

end Compact

end OrbitModel

/-! ## 4. 具体層(有限次) —— mathlib の正規底からの witness

★ここが**水準ごとの条件の非空虚性**である:有限次 Galois 拡大 `E/F` には、
`Φ_λ` が全単射になる汎関数 `λ` が実際に存在する(正規底の第 1 座標)。 -/

section FiniteGalois

open OrbitModel

variable (F E : Type*) [Field F] [Field E] [Algebra F E] [FiniteDimensional F E] [IsGalois F E]

/-- 正規底の「第 1 座標」—— **双対正規底生成元**。 -/
noncomputable def normalBasisCoord : E →+ F :=
  ((IsGalois.normalBasis F E).coord 1).toAddMonoidHom

theorem orbitHom_normalBasisCoord (x : E) (σ : E ≃ₐ[F] E) :
    orbitHom (G := E ≃ₐ[F] E) (normalBasisCoord F E) x σ
      = (IsGalois.normalBasis F E).equivFun x σ := by
  classical
  show (IsGalois.normalBasis F E).coord 1 (σ⁻¹ • x) = _
  rw [Module.Basis.coord_apply, ← Module.Basis.equivFun_apply]
  have h := normalBasis_equivFun_smul F E σ⁻¹ x 1
  rw [show ((σ⁻¹ : E ≃ₐ[F] E)⁻¹ * 1) = σ by simp] at h
  exact h

/-- ★★**有限次 Galois の水準では仮説が満たされる**(非空虚性、`sorry` 無し)。

`Φ_λ` は正規底での座標表示そのものになる。 -/
theorem orbitHom_normalBasisCoord_bijective :
    Function.Bijective (orbitHom (G := E ≃ₐ[F] E) (normalBasisCoord F E)) := by
  classical
  have h : (orbitHom (G := E ≃ₐ[F] E) (normalBasisCoord F E) : E → ((E ≃ₐ[F] E) → F))
      = (IsGalois.normalBasis F E).equivFun := by
    funext x
    funext σ
    exact orbitHom_normalBasisCoord F E x σ
  rw [h]
  exact (IsGalois.normalBasis F E).equivFun.bijective

end FiniteGalois


/-! ## 5. 具体層(無限次) —— `Γ_K` の水準への分解 -/

section Concrete

open ABC3.Skeleton.PGC OrbitModel

variable {p : ℕ} [Fact p.Prime]

/-- `K.carrier` は標数 0(`ℚ_p` 上の代数)。★以下では局所インスタンスとして使う。 -/
theorem carrier_charZero (K : PAdicLocalField p) : CharZero K.carrier :=
  charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective

attribute [local instance] carrier_charZero

variable (K : PAdicLocalField p)

/-- 水準 `L` の包含 `L →+ K̄`。 -/
def levelIota (L : IntermediateField K.carrier K.closure) : L →+ K.closure where
  toFun b := (b : K.closure)
  map_zero' := rfl
  map_add' _ _ := rfl

@[simp] theorem levelIota_apply (L : IntermediateField K.carrier K.closure) (b : L) :
    levelIota K L b = (b : K.closure) := rfl

theorem levelIota_injective (L : IntermediateField K.carrier K.closure) :
    Function.Injective (levelIota K L) := fun _ _ h => Subtype.ext h

/-- 水準の射影 `Γ_K ↠ Gal(L/K)`。 -/
noncomputable def levelPi (L : IntermediateField K.carrier K.closure) [Normal K.carrier L] :
    K.absGal →* (L ≃ₐ[K.carrier] L) := AlgEquiv.restrictNormalHom (L : Type _)

theorem coe_levelPi_apply (L : IntermediateField K.carrier K.closure) [Normal K.carrier L]
    (g : K.absGal) (b : L) : ((levelPi K L g b : L) : K.closure) = g (b : K.closure) := by
  have hc := AlgEquiv.restrictNormal_commutes g (L : Type _) b
  simp only [levelPi, AlgEquiv.restrictNormalHom, MonoidHom.mk'_apply]
  exact hc

theorem levelPi_surjective (L : IntermediateField K.carrier K.closure) [Normal K.carrier L] :
    Function.Surjective (levelPi K L) := AlgEquiv.restrictNormalHom_surjective K.closure

/-- ★遷移の同変性(具体層)。 -/
theorem levelIota_smul (L : IntermediateField K.carrier K.closure) [Normal K.carrier L]
    (g : K.absGal) (b : L) :
    levelIota K L (levelPi K L g • b) = g • levelIota K L b :=
  coe_levelPi_apply K L g b

/-- `ker (Γ_K ↠ Gal(L/K)) = L.fixingSubgroup`。 -/
theorem levelPi_eq_one_iff (L : IntermediateField K.carrier K.closure) [Normal K.carrier L]
    (g : K.absGal) : levelPi K L g = 1 ↔ g ∈ L.fixingSubgroup := by
  rw [IntermediateField.mem_fixingSubgroup_iff]
  constructor
  · intro h x hx
    have hc := coe_levelPi_apply K L g ⟨x, hx⟩
    rw [h] at hc
    exact hc.symm
  · intro h
    ext b
    have hc := coe_levelPi_apply K L g b
    rw [h (b : K.closure) b.2] at hc
    exact hc

/-! ### 5.2 被覆(★水準への分解に必要な 2 つの事実) -/

/-- ★(a) 有限個の元は 1 つの有限次 Galois 部分拡大に入る。 -/
theorem exists_finiteGalois_of_finite (S : Set K.closure) (hS : S.Finite) :
    ∃ L : IntermediateField K.carrier K.closure, FiniteDimensional K.carrier L ∧
      IsGalois K.carrier L ∧ S ⊆ (L : Set K.closure) := by
  haveI : Finite S := hS.to_subtype
  exact ⟨(FiniteGaloisIntermediateField.adjoin K.carrier S).toIntermediateField,
    inferInstance, inferInstance, FiniteGaloisIntermediateField.subset_adjoin K.carrier S⟩

/-- ★(b) `Γ_K` 上の局所定数関数は有限次 Galois の水準を経由する。

★抽象核 III(コンパクト群の一様性)+ krull 位相の 1 の近傍基底、の 2 段。 -/
theorem exists_level_factor {M : Type*} (f : LocallyConstant K.absGal M) :
    ∃ (L : IntermediateField K.carrier K.closure) (_ : FiniteDimensional K.carrier L)
      (_ : IsGalois K.carrier L) (f₀ : (L ≃ₐ[K.carrier] L) → M),
      ∀ g : K.absGal, f g = f₀ (levelPi K L g) := by
  classical
  obtain ⟨U, hU, hUspec⟩ := OrbitModel.exists_nhds_one_forall_mul f
  obtain ⟨L, hL⟩ := (InfiniteGalois.krullTopology_mem_nhds_one_iff_of_isGalois U).1 hU
  refine ⟨L.toIntermediateField, inferInstance, inferInstance,
    fun σ => f (Function.surjInv (levelPi_surjective K L.toIntermediateField) σ), ?_⟩
  intro g
  have hpi : levelPi K L.toIntermediateField
      (Function.surjInv (levelPi_surjective K L.toIntermediateField) (levelPi K _ g))
      = levelPi K L.toIntermediateField g :=
    Function.surjInv_eq (levelPi_surjective K L.toIntermediateField) _
  have hker : levelPi K L.toIntermediateField
      (g⁻¹ * Function.surjInv (levelPi_surjective K L.toIntermediateField) (levelPi K _ g)) = 1 := by
    rw [map_mul, map_inv, hpi]
    simp
  have hmem := hL ((levelPi_eq_one_iff K L.toIntermediateField _).1 hker)
  have hval := hUspec g _ hmem
  rw [mul_inv_cancel_left] at hval
  exact hval.symm

/-! ### 5.3 ★★★還元 -/

/-- ★★★**水準ごとの条件から `SmoothModelCarrier K`**。

「1 つの汎関数 `λ : K̄ →+ K` であって、**すべての**有限次 Galois 部分拡大 `L/K` について
`Φ_{λ|_L}` が全単射」が取れれば、`K̄ ≅ C^∞(Γ_K, K)` が `Γ_K`-同変に従う。

★塔をまたぐ両立性(coherence)は**要求されない** —— `λ` を制限するだけだから。 -/
theorem smoothModelCarrier_of_levelwise (lam : K.closure →+ K.carrier)
    (hlev : ∀ (L : IntermediateField K.carrier K.closure), FiniteDimensional K.carrier L →
      IsGalois K.carrier L →
      Function.Bijective (orbitHom (G := L ≃ₐ[K.carrier] L) (lam.comp (levelIota K L)))) :
    SmoothModelCarrier K := by
  have hlc : ∀ a : K.closure, IsLocallyConstant (orbitFun (G := K.absGal) lam a) :=
    fun a => isLocallyConstant_orbitFun lam a (stabilizer_isOpen_of_isIntegral a)
  refine ⟨semilinearOfBijective lam hlc ⟨?_, ?_⟩⟩
  · intro a₁ a₂ h
    obtain ⟨L, hfd, hgal, hsub⟩ := exists_finiteGalois_of_finite K {a₁, a₂} (Set.toFinite _)
    haveI := hfd
    haveI := hgal
    have h1 : a₁ ∈ L := hsub (by simp)
    have h2 : a₂ ∈ L := hsub (by simp)
    have hfun : orbitFun (G := K.absGal) lam (levelIota K L ⟨a₁, h1⟩)
        = orbitFun (G := K.absGal) lam (levelIota K L ⟨a₂, h2⟩) := by
      funext g
      exact congrFun (congrArg (fun u : LocallyConstant K.absGal K.carrier =>
        (u : K.absGal → K.carrier)) h) g
    have hb := eq_of_orbitFun_eq_of_level (levelPi K L) (levelPi_surjective K L) (levelIota K L)
      (levelIota_smul K L) lam (hlev L hfd hgal).1 ⟨a₁, h1⟩ ⟨a₂, h2⟩ hfun
    exact congrArg Subtype.val hb
  · intro f
    obtain ⟨L, hfd, hgal, f₀, hf⟩ := exists_level_factor K f
    haveI := hfd
    haveI := hgal
    obtain ⟨b, hb⟩ := exists_orbitFun_eq_of_level (levelPi K L) (levelIota K L)
      (levelIota_smul K L) lam (hlev L hfd hgal).2 f₀
    refine ⟨levelIota K L b, ?_⟩
    ext g
    show orbitFun (G := K.absGal) lam (levelIota K L b) g = f g
    rw [hb]
    exact (hf g).symm

/-! ### 5.4 ★残っている唯一の穴に名前を付ける -/

/-- ★★**両立汎関数の存在** —— `SmoothModelCarrier K` に残る唯一の内容。

「**1 つの** `λ : K̄ →+ K` であって、**すべての**有限次 Galois 部分拡大 `L/K` について
`Φ_{λ|_L} : L → Maps(Gal(L/K), K)` が全単射」。

★各水準**単独**では必ず満たせる(`exists_levelwise_functional`、mathlib の正規底定理)。
難しいのは「1 つの `λ` で全水準を同時に」であり、そこに
(i) Maschke による延長(`OrbitModel.exists_unit_of_central_idempotent` が受け皿)と
(ii) 有限次部分拡大の可算性(Krasner)が要る。 -/
def HasCoherentFunctional (K : PAdicLocalField p) : Prop :=
  ∃ lam : K.closure →+ K.carrier, ∀ (L : IntermediateField K.carrier K.closure),
    FiniteDimensional K.carrier L → IsGalois K.carrier L →
    Function.Bijective (orbitHom (G := L ≃ₐ[K.carrier] L) (lam.comp (levelIota K L)))

theorem smoothModelCarrier_of_hasCoherentFunctional (K : PAdicLocalField p)
    (h : HasCoherentFunctional K) : SmoothModelCarrier K := by
  obtain ⟨lam, hlam⟩ := h
  exact smoothModelCarrier_of_levelwise K lam hlam

/-- ★★★**[pGC] Proposition 2.1 への最終配線**。

`HasCoherentFunctional` が全ての `K` で成り立てば、`K̄` は `Γ_K` から群論的に回復される。 -/
theorem prop_2_1_of_hasCoherentFunctional
    (h : ∀ K : PAdicLocalField p, HasCoherentFunctional K) :
    RecoverableAsAddModule (p := p) (fun K => K.closure) :=
  prop_2_1_of_smoothModelCarrier (fun K => smoothModelCarrier_of_hasCoherentFunctional K (h K))

/-- ★**非空虚性(水準ごと)**:各有限次 Galois 部分拡大には良い汎関数が実際に存在する
(mathlib の正規底定理)。★ただしここで得られる `μ` は水準ごとに別物であり、
`HasCoherentFunctional` が要求する「1 つの `λ` の制限」にはなっていない —— そこが穴である。 -/
theorem exists_levelwise_functional (L : IntermediateField K.carrier K.closure)
    [FiniteDimensional K.carrier L] [IsGalois K.carrier L] :
    ∃ mu : L →+ K.carrier, Function.Bijective (orbitHom (G := L ≃ₐ[K.carrier] L) mu) :=
  ⟨normalBasisCoord K.carrier L, orbitHom_normalBasisCoord_bijective K.carrier L⟩

/-! ### 5.5 ★★逆向き —— 還元に**損失が無い**こと

`SmoothModelCarrier K` から `HasCoherentFunctional K` が戻る。すなわち
§5.3 で切り出した仮説は、目標と**同値**である(強すぎも弱すぎもしない)。

★鍵は Frobenius 相互律:`Γ_K`-同変な加法同型 `K̄ ≅ C^∞(Γ_K, K)` は
**必ず** `Φ_λ` の形をしている(`λ := (1 での評価) ∘ E`)。 -/

/-- `ker (Γ_K ↠ Gal(L/K))` は開(= `L.fixingSubgroup`)。 -/
theorem isOpen_ker_levelPi (L : IntermediateField K.carrier K.closure)
    [FiniteDimensional K.carrier L] [Normal K.carrier L] :
    IsOpen (((levelPi K L).ker : Subgroup K.absGal) : Set K.absGal) := by
  have h : (levelPi K L).ker = L.fixingSubgroup := by
    ext g
    rw [MonoidHom.mem_ker]
    exact levelPi_eq_one_iff K L g
  rw [h]
  exact L.fixingSubgroup_isOpen

/-- 水準の核で固定される元は水準に属する(Galois 降下)。 -/
theorem mem_range_levelIota_of_fixed (L : IntermediateField K.carrier K.closure)
    [FiniteDimensional K.carrier L] [IsGalois K.carrier L] (a : K.closure)
    (h : ∀ n : K.absGal, levelPi K L n = 1 → n • a = a) : a ∈ Set.range (levelIota K L) := by
  have ha : a ∈ IntermediateField.fixedField L.fixingSubgroup := by
    rw [IntermediateField.mem_fixedField_iff]
    intro g hg
    exact h g ((levelPi_eq_one_iff K L g).2 hg)
  rw [InfiniteGalois.fixedField_fixingSubgroup] at ha
  exact ⟨⟨a, ha⟩, rfl⟩

/-- ★★★**逆向きの含意**。 -/
theorem hasCoherentFunctional_of_smoothModelCarrier (K : PAdicLocalField p)
    (h : SmoothModelCarrier K) : HasCoherentFunctional K := by
  obtain ⟨E⟩ := h
  refine ⟨{ toFun := fun a => E.toAddEquiv a 1
            map_zero' := by simp
            map_add' := fun x y => by simp }, ?_⟩
  set lam : K.closure →+ K.carrier :=
    { toFun := fun a => E.toAddEquiv a 1
      map_zero' := by simp
      map_add' := fun x y => by simp } with hlam
  have hEeq : ∀ (a : K.closure) (g : K.absGal), lam (g⁻¹ • a) = E.toAddEquiv a g := by
    intro a g
    have hs := E.map_smul g⁻¹ a
    show E.toAddEquiv (g⁻¹ • a) 1 = E.toAddEquiv a g
    rw [hs]
    show E.toAddEquiv a ((g⁻¹)⁻¹ * 1) = E.toAddEquiv a g
    rw [inv_inv, mul_one]
  have hlc : ∀ a : K.closure, IsLocallyConstant (orbitFun (G := K.absGal) lam a) :=
    fun a => isLocallyConstant_orbitFun lam a (stabilizer_isOpen_of_isIntegral a)
  have hbij : Function.Bijective (orbitLC lam hlc) := by
    have heq : (orbitLC lam hlc : K.closure → LocallyConstant K.absGal K.carrier)
        = E.toAddEquiv := by
      funext a
      ext g
      exact hEeq a g
    rw [heq]
    exact E.toAddEquiv.bijective
  intro L hfd hgal
  haveI := hfd
  haveI := hgal
  exact bijective_descent_lc (levelPi K L) (levelPi_surjective K L) (levelIota K L)
    (levelIota_injective K L) (levelIota_smul K L) (mem_range_levelIota_of_fixed K L)
    (fun f => isLocallyConstant_comp_of_isOpen_ker (levelPi K L) (isOpen_ker_levelPi K L) f)
    lam hlc hbij

/-- ★★★**同値**:`SmoothModelCarrier K ↔ HasCoherentFunctional K`。 -/
theorem smoothModelCarrier_iff_hasCoherentFunctional (K : PAdicLocalField p) :
    SmoothModelCarrier K ↔ HasCoherentFunctional K :=
  ⟨hasCoherentFunctional_of_smoothModelCarrier K, smoothModelCarrier_of_hasCoherentFunctional K⟩

end Concrete

end ABC3.Found.PGC
