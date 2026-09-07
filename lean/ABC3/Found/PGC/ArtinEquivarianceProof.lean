import ABC3.Found.PGC.ArtinEquivariance
import ABC3.Found.PGC.TwistedLTComplete

/-!
# 経路 Λ13 —— 残った壁を「Artin 写像の**単数部分**の `σ`-同変性」に落とす

[pGC] Proposition 1.1 への経路 Λ の最後の 1 本(Λ12 `ArtinEquivariance.lean`)は、壁を

  「`Λ_n(S) ≅ μ_{p^n}` を `Γ_F`-同変に取れるか」

という 5 つの同値な言い方に整えた。本ファイルはその壁を**さらに 1 段**落として、

  ★**`ArtinUnitEquivariance`** ——「`Gal(L^{ab}/L) ≅ 𝒪_L^× × Ẑ`(局所類体論の分解)の
  **第 1 成分**が `σ_g` と同変である」

にする。★★**そしてそれが元の壁と同値であることを証明する**
(`artinUnitEquivariance_iff_cyclotomeConj`)。

## ★★何が新しく分かったか(★測定結果。楽観的に書かない)

| 段 | 内容 | 本ファイルの結果 |
|---|---|---|
| (a) | 壁を Λ9 の座標(`𝒪_L^× × Ẑ`)に翻訳する | ★**証明した**(`cyclotomeConjIsCyclotomic_of_artinUnitEquivariance`) |
| (b) | その翻訳が**同値**であること | ★**証明した**(`artinUnitEquivariance_iff_cyclotomeConj`) |
| (c) | `g ∈ S` の場合(壁の易しい半分) | ★**証明した**(`exists_artinUnitEquivariance_of_mem`) |
| (d) | `σ_g(π)` が別の素元であること | ★**証明した**(`maximalIdeal_eq_span_integerRingEquiv`) |
| (e) | `f^{σ_g}` が `σ_g(π)` の Lubin-Tate 級数であること | ★**証明した**(`lubinTateSeries_integerRingEquiv`) |
| (f) | 「点で評価する層」(`σ` が `Λ_{f,n}` を `Λ_{f^σ,n}` に運ぶ) | ★**証明した**(`map_mem_iteratedLubinTateTorsionPoints`) |
| (g) | `g ∉ S` の場合の同変性そのもの | ★**残った**(`ArtinUnitEquivariance`) |

★★**同変性は出ていない。** 出たのは「壁の位置を 1 段下げたこと」と、
その 1 段下で使う道具 (d)(e)(f) である。

## ★★残った壁は何か(★次の節点への申し送り、名指しで)

```
ArtinUnitEquivariance p :
  ∀ F n S (hS : S.Normal) (hopen : IsOpen S), S ≤ muFixer F (p^n) → ∀ g : Γ_F,
    ∃ E : Gal(L^{ab}/L) ≃* 𝒪_L^× × Ẑ,
      ∀ x, x^{p^n} = 1 → u(E(g x g^{-1})) = σ_g (u(E x))
```
ここで `L := L_S`、`σ_g := fixedFieldAut S g`、`u` は `𝒪_L^× → L^×` である。

★★**これは局所類体論の「Artin 写像の関手性」そのもの**であって、それ以上でも以下でもない。
★★**Ẑ 成分(Frobenius 部分)は結論に現れない** —— 壁は `𝒪_L^×` 成分だけである。
これは本ファイルで**測って分かったこと**である
(Λ12 の申し送りは「Artin 写像の同変性」としか書いていなかった)。

★Lubin-Tate で埋めるなら、残っているのは

1. `σ_g` は `𝒪_L` の環自己同型を誘導する —— ★**本ファイルで済ませた**
   (`fixedFieldIntegerAut`、在庫 `integerRingEquiv` に代入するだけ)。
2. `σ_g(π)` は別の素元、`f^{σ_g}` はその Lubin-Tate 級数 —— ★**本ファイルで済ませた**
   ((d)(e))。したがって `K_{σ_g(π)}` を作る材料は揃っている。
3. `σ_g` は `Λ_{f,m}` を `Λ_{f^{σ_g},m}` に写す(★**点で評価する層**) ——
   ★★**本ファイルで済ませた**(`map_mem_iteratedLubinTateTorsionPoints`、§13)。
   経路は「§3 の抽象核 C(半線型な `Φ` は係数をひねった多項式の値を運ぶ)+
   §11 の抽象核 D(Weierstrass 標準分解の自然性)+ §12 の `D_n` の自然性」。
4. ★★**残っているのはここだけ**:`σ_g` が `K_{f}` を `K_{f^{σ_g}}` に写すこと
   (§13 を `adjoin` に持ち上げる段)と、
   `ρ_{f^{σ_g}}(g τ g^{-1}) = σ_g(ρ_f(τ))`、そして素元非依存性
   (`reciprocityHom_eq_of_intertwiner` /
   `lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer`)で `π` に戻す段。
   ★**これが次の節点である。**

## ★★退化の自己検査

1. **空虚ではない**。`natCard_cyclotome_eq`(Λ11)と `natCard_rootsOfUnity_closure`(Λ12)が
   両辺の位数をちょうど `p^n` と与えている。`ArtinUnitEquivariance` は
   `CyclotomeConjIsCyclotomic` と**同値**(`artinUnitEquivariance_iff_cyclotomeConj`)なので、
   「通すために弱めた」形ではない。
2. **易しい半分は実際に成り立っている**。`g ∈ S` なら
   `abelianGalConj_of_mem`(左辺が恒等)と `fixedFieldAut_of_mem`(右辺が恒等)が
   **独立の理由で**両辺を恒等にする(`exists_artinUnitEquivariance_of_mem`)。
   ★したがって壁の内容は `Γ_F ⧸ S = Gal(L_S/F)` の上にある。
3. **`S ≤ muFixer F (p^n)` を落としていない**。落とすと同値性の逆向きが壊れる
   ——`σ_g` が `μ_{p^n}(L_S)` の上で `χ_n(g)` 倍であること
   (`rootsOfUnityFixedFieldConj_eq_pow`)がこの仮定を使う。
4. **式の向きを変えていない**。Λ9 の `hΨ`(共変、`Art(σ x) = σ Art(x) σ^{-1}`)と
   同じ向きである。★Λ8 の符号の測定(§2.2 の Weil 群は**幾何** Frobenius の巾)は
   `Ẑ` 成分の話であり、本ファイルの壁には現れない —— したがって壊していない。
5. **`n = 0` でも成り立つ**(両辺とも自明群)。

## 逸脱の記録

* 新しい仮定を 1 つも置いていない。`ArtinUnitEquivariance` は
  `CyclotomeConjIsCyclotomic` と**同値**である。
* 原典(pGC §1)の論拠は局所 Tate 双対性であって局所類体論ではない。この逸脱は
  Λ10/Λ11/Λ12 と同じく `ResearchPaper/pgc-goal.md` に記録済み。
* (d)(e) は `𝒪_L` の自己同型を `integerRingEquiv`(`ResidueCardinality.lean`、
  ノルム保存から作られている)経由で取っている。★`σ_g` が `𝒪_L` を保つことを
  「付値を保つ」ではなく「スペクトルノルムを保つ」で出しているが、
  局所体では同じことである(在庫 `norm_algEquiv_carrier`)。

## ★測定の記録(★次の agent が再実行できる形で)

* ★**mathlib を先に引いた**。`grep -n "IsLocalHom" .cache/mathlib-index.txt | grep -iE "equiv"`
  で `isLocalHom_equiv`(`Algebra/Group/Units/Equiv.lean:214`、`@[instance]`)が出た
  ——「局所環の同型は局所準同型」を自分で書かずに済んだ。
  `grep -n "ResidueField.map" .cache/mathlib-index.txt` で
  `IsLocalRing.ResidueField.map_comp_residue` / `mapEquiv` が出た。
  `grep -n "weierstrassDistinguished" .cache/mathlib-index.txt` で
  `PowerSeries.IsWeierstrassFactorization.unique`(分解の一意性)が出た ——
  ★**§11 の自然性はこれ 1 本で出る**(自分で一意性を書いていない)。
  ★**「mathlib に無い」と書いた箇所は本ファイルに 1 つも無い。**
* ★**#158(同名で書いて `already been declared`)**: `abelianGalConj` /
  `fixedFieldIntegerAut` / `hom_eval₂_twist` の 3 つを試して 0 件。
  一方 `grep -rn "integerRingEquiv" lean/ABC3/` は
  `integerRingEquiv`(`ResidueCardinality.lean:77`)を**引き当てた** ——
  ★「`σ` が整数環を保つ」段は**既に在庫にあった**ので書かずに済んだ。
  ★同じく `grep -rn "map_subst\|subst_map" lean/ABC3/Found/PGC/*.lean` が
  `map_iteratedLubinTate`(`DworkThetaEval.lean:278`、係数写像は自己合成と可換)と
  `map_subst_powerSeries`(`DworkThetaStep2.lean:258`)を**引き当てた** ——
  ★§12 はこの 2 本のおかげで 15 行で済んだ。★★**「作られた目的」ではなく
  「型」で引くと当たる**の実例がまた 2 件増えた(どちらも Dwork の θ のために
  作られたもので、捩れ点のためではない)。
* ★**素元非依存性(`LubinTateUniformizerIndependence`)は本ファイルでも消費しなかった。**
  消費するのは上の 4 である。
* ★**`LubinTateEndoTwisted`(冪級数のねじれ版)も消費しなかった。**
  本ファイルが作った (d)(e) はその**入力**(`f^{σ}` が Lubin-Tate 級数であること)であり、
  `[θ]_{f,f^σ}` を実際に使うのは上の 4 である。
* ★**#59 には当たっていない。** 中間体は `L_S = fixedField S` の 1 層しか出てこない
  (定型 (b):`Γ_F` から 1 層で書く)。`abelianClosure L` は `L` の中で閉じている。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

open scoped NormedField Valued Classical IsMulCommutative

/-! ## 1. ★抽象核 A —— 同型に沿った「冪作用」の輸送(純群論)

★この節に体・付値・分岐・位相の語彙は 1 つも出てこない。 -/

section AbstractCore

variable {A B : Type*} [Group A] [Group B]

/-- **★抽象核 A1** —— 同型 `α` が `Ψ` を `Ψ'` に移し、`Ψ'` が `c` 乗なら、`Ψ` も `c` 乗。

★`α` の単射性しか使わない。 -/
theorem eq_pow_of_transport (α : A ≃* B) (Ψ : A ≃* A) (Ψ' : B ≃* B) (c : ℕ)
    (hcomm : ∀ x, α (Ψ x) = Ψ' (α x)) (h : ∀ y, Ψ' y = y ^ c) (x : A) :
    Ψ x = x ^ c := by
  apply α.injective
  rw [hcomm, h, ← map_pow]

/-- **★抽象核 A1'(逆向き)** —— `Ψ` が `c` 乗なら `Ψ'` も `c` 乗。

★`α` の全射性しか使わない。★A1 と合わせて「輸送は情報を落とさない」。 -/
theorem transport_eq_pow (α : A ≃* B) (Ψ : A ≃* A) (Ψ' : B ≃* B) (c : ℕ)
    (hcomm : ∀ x, α (Ψ x) = Ψ' (α x)) (h : ∀ x, Ψ x = x ^ c) (y : B) :
    Ψ' y = y ^ c := by
  obtain ⟨x, rfl⟩ := α.surjective y
  rw [← hcomm, h, map_pow]

end AbstractCore

/-! ## 2. ★抽象核 B —— 局所環の自己同型は Lubin-Tate 級数を Lubin-Tate 級数に移す

★この節に体・Galois・分岐の語彙は 1 つも出てこない。一般の可換局所環である。
★★これが Λ12 の申し送り (a)「`σ_g(π)` が別の素元であること」の抽象核である。 -/

section TwistCore

variable {A : Type*} [CommRing A] [IsLocalRing A]

/-- **★抽象核 B1** —— 局所環の自己同型は素元を素元に移す。

`𝔪 = (π)` なら `𝔪 = (φ π)`。★極大イデアルは「単元でない元の集合」という
**環同型で不変な**記述を持つので、これは付値を 1 度も使わずに出る。 -/
theorem maximalIdeal_eq_span_map (φ : A ≃+* A) {π : A}
    (h : IsLocalRing.maximalIdeal A = Ideal.span {π}) :
    IsLocalRing.maximalIdeal A = Ideal.span {φ π} := by
  ext x
  rw [Ideal.mem_span_singleton, IsLocalRing.mem_maximalIdeal]
  constructor
  · intro hx
    have hx' : φ.symm x ∈ IsLocalRing.maximalIdeal A := by
      rw [IsLocalRing.mem_maximalIdeal]
      intro hu
      exact hx (by simpa using hu.map φ)
    rw [h, Ideal.mem_span_singleton] at hx'
    have hdvd := map_dvd (φ : A →+* A) hx'
    simpa using hdvd
  · rintro ⟨c, rfl⟩
    intro hu
    have hπ : π ∈ IsLocalRing.maximalIdeal A := by
      rw [h, Ideal.mem_span_singleton]
    have hφπ : IsUnit (φ π) := isUnit_of_mul_isUnit_left hu
    exact (IsLocalRing.mem_maximalIdeal π).mp hπ (by simpa using hφπ.map φ.symm)

omit [IsLocalRing A] in
/-- ★抽象核 B2(定数項) —— 係数をひねっても `f(0) = 0` は保たれる。 -/
theorem coeff_zero_map_twist (φ : A →+* A) {f : PowerSeries A}
    (hf0 : PowerSeries.coeff 0 f = 0) : PowerSeries.coeff 0 (PowerSeries.map φ f) = 0 := by
  rw [PowerSeries.coeff_map, hf0, map_zero]

omit [IsLocalRing A] in
/-- ★抽象核 B3(1 次の係数) —— `f'(0) = π` は `f'(0) = φ π` になる。 -/
theorem coeff_one_map_twist (φ : A →+* A) {f : PowerSeries A} {π : A}
    (hf1 : PowerSeries.coeff 1 f = π) : PowerSeries.coeff 1 (PowerSeries.map φ f) = φ π := by
  rw [PowerSeries.coeff_map, hf1]

/-- ★★抽象核 B4(剰余体の条件) —— `f ≡ X^q (mod 𝔪)` は係数をひねっても保たれる。

`φ` は剰余体の自己同型 `φ̄` を誘導し、`residue ∘ φ = φ̄ ∘ residue`
(mathlib `IsLocalRing.ResidueField.map_comp_residue`)。`φ̄` は `X^q` を `X^q` に写す。 -/
theorem map_residue_map_twist (φ : A ≃+* A) {f : PowerSeries A} {q : ℕ}
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ q) :
    PowerSeries.map (IsLocalRing.residue A) (PowerSeries.map (φ : A →+* A) f)
      = PowerSeries.X ^ q := by
  calc PowerSeries.map (IsLocalRing.residue A) (PowerSeries.map (φ : A →+* A) f)
      = PowerSeries.map ((IsLocalRing.residue A).comp (φ : A →+* A)) f := by
        rw [PowerSeries.map_comp]; rfl
    _ = PowerSeries.map ((IsLocalRing.ResidueField.map (φ : A →+* A)).comp
          (IsLocalRing.residue A)) f := by
        rw [IsLocalRing.ResidueField.map_comp_residue]
    _ = PowerSeries.map (IsLocalRing.ResidueField.map (φ : A →+* A))
          (PowerSeries.map (IsLocalRing.residue A) f) := by
        rw [PowerSeries.map_comp]; rfl
    _ = PowerSeries.X ^ q := by rw [hf, map_pow, PowerSeries.map_X]

end TwistCore

/-! ## 3. ★抽象核 C —— 「点で評価する層」

★★Λ12 の申し送り (b)「`σ` が捩れ点を運ぶ層」の**抽象核**である。
★この節にも体・付値・分岐・Lubin-Tate の語彙は 1 つも出てこない。 -/

section EvalTwist

variable {A B : Type*} [CommRing A] [CommRing B]

/-- **★★抽象核 C1(点で評価する層)** ——
`Φ` が `φ` について半線型(`Φ ∘ ι = ι ∘ φ`)なら、`Φ` は
「係数を `φ` でひねった多項式」の値を運ぶ:

  `Φ (P(x)) = (P^φ)(Φ x)`。

★これが「`σ` は `f`-捩れ点を `f^σ`-捩れ点に写す」の**核**である。 -/
theorem hom_eval₂_twist (φ : A →+* A) (Φ : B →+* B) (ι : A →+* B)
    (h : ∀ a, Φ (ι a) = ι (φ a)) (P : Polynomial A) (x : B) :
    Φ (Polynomial.eval₂ ι x P) = Polynomial.eval₂ ι (Φ x) (Polynomial.map φ P) := by
  rw [Polynomial.hom_eval₂, Polynomial.eval₂_map]
  congr 1
  exact RingHom.ext h

/-- ★★抽象核 C2 —— 根は「ひねった多項式」の根に写る。 -/
theorem eval₂_map_twist_eq_zero (φ : A →+* A) (Φ : B →+* B) (ι : A →+* B)
    (h : ∀ a, Φ (ι a) = ι (φ a)) (P : Polynomial A) {x : B}
    (hx : Polynomial.eval₂ ι x P = 0) :
    Polynomial.eval₂ ι (Φ x) (Polynomial.map φ P) = 0 := by
  rw [← hom_eval₂_twist φ Φ ι h P x, hx, map_zero]

end EvalTwist

variable {p : ℕ} [Fact p.Prime]

/-! ## 4. 具体層 1 —— `σ` が `𝒪_K` に誘導する自己同型と、ひねった Lubin-Tate 級数

★★Λ12 の申し送り (a) の具体層。★在庫 `integerRingEquiv`
(`ResidueCardinality.lean`)を §2 の抽象核に代入するだけである。 -/

section Uniformizer

variable (K : PAdicLocalField p) (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier)

/-- **★`σ(π)` は別の素元** —— `𝔪 = (π)` なら `𝔪 = (σ π)`。

★★Λ12 が名指しした入力 (a) の前半。★`σ` はノルムを保つので `𝒪_K` を保ち
(`integerRingEquiv`)、その自己同型は極大イデアルを保つ(§2 の抽象核 B1)。 -/
theorem maximalIdeal_eq_span_integerRingEquiv {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) :
    IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {integerRingEquiv σ π} :=
  maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax

/-- ★`σ(π) ≠ 0`。 -/
theorem integerRingEquiv_ne_zero {π : 𝒪[K.carrier]} (hπne0 : π ≠ 0) :
    integerRingEquiv σ π ≠ 0 := by
  intro h
  exact hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero]))

/-- **★★`f^σ` は `σ(π)` の Lubin-Tate 級数** —— Λ12 が名指しした入力 (a) の後半。

★★これで「`σ` の像側の Lubin-Tate 塔 `K_{σ(π)}`」を作る材料が揃う
(`lubinTateClosure K hq (maximalIdeal_eq_span_integerRingEquiv …) … (PowerSeries.map … f) …`)。
★**塔が一致すること**(素元非依存性)はここでは主張していない —— それは
`lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer` の仕事である。 -/
theorem lubinTateSeries_integerRingEquiv
    {π : 𝒪[K.carrier]} {q : ℕ} (f : PowerSeries 𝒪[K.carrier])
    (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ q) :
    PowerSeries.coeff 0
        (PowerSeries.map (integerRingEquiv σ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f) = 0
      ∧ PowerSeries.coeff 1
        (PowerSeries.map (integerRingEquiv σ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          = integerRingEquiv σ π
      ∧ PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier])
        (PowerSeries.map (integerRingEquiv σ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          = PowerSeries.X ^ q :=
  ⟨coeff_zero_map_twist _ hf0, coeff_one_map_twist _ hf1,
    map_residue_map_twist (integerRingEquiv σ) hf⟩

end Uniformizer

/-- **`σ_g` が誘導する `𝒪_{L_S}` の環自己同型**。

★`fixedFieldAut S g` は `F.carrier` 上の代数同型なので、`ℚ_p` 上に制限してから
`integerRingEquiv` に代入する。★中間体は `L_S` の 1 層しか出てこない(#59 の定型 (b))。 -/
noncomputable def fixedFieldIntegerAut (F : PAdicLocalField p) (S : Subgroup F.absGal)
    [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal) :
    𝒪[(fixedFieldLocalField F S hopen).carrier] ≃+* 𝒪[(fixedFieldLocalField F S hopen).carrier] :=
  integerRingEquiv ((fixedFieldAut S g).restrictScalars ℚ_[p])

@[simp] theorem coe_fixedFieldIntegerAut (F : PAdicLocalField p) (S : Subgroup F.absGal)
    [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (z : 𝒪[(fixedFieldLocalField F S hopen).carrier]) :
    ((fixedFieldIntegerAut F S hopen g z : 𝒪[(fixedFieldLocalField F S hopen).carrier])
        : (fixedFieldLocalField F S hopen).carrier)
      = fixedFieldAut S g (z : (fixedFieldLocalField F S hopen).carrier) := rfl

/-! ## 5. 具体層 2 —— `Gal(L^{ab}/L)` の上の `σ_g`-半線型共役 `Ψ_g` -/

/-- **★`Gal(L_S^{ab}/L_S)` の上の共役作用 `Ψ_g`**。

Λ12 の `absGalConjCME`(`Γ_{L_S}` の上の `τ ↦ g τ g^{-1}`)を、Λ11 の
`topAbelianizationEquivAbelianGal`(`Γ_L^{ab} ≅ Gal(L^{ab}/L)`)で移したもの。

★★**Artin 写像の同変性 `Art_L(σ_g x) = g Art_L(x) g^{-1}` の右辺に現れる作用がこれ**である。 -/
noncomputable def abelianGalConj (F : PAdicLocalField p) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal) :
    Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
      ≃*
    Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier) :=
  ((topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)).symm.trans
      (topAbelianizationCME (absGalConjCME F S hopen g)).toMulEquiv).trans
    (topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen))

theorem abelianGalConj_apply (F : PAdicLocalField p) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (x : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)) :
    abelianGalConj F S hopen g x
      = topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)
          (topAbelianizationCME (absGalConjCME F S hopen g)
            ((topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)).symm x)) := rfl

/-- ★**`Ψ_g` は `Γ_{L_S}` の側の共役と対応する**(捩れの上で)。

★`cyclotome A m` と `powTorsion (A^{ab}) m` は**定義が同じ**なので、
橋は `topAbelianizationEquivAbelianGal` の 1 本で足りる。 -/
theorem powTorsionCongr_abelianGalConj (F : PAdicLocalField p) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) (m : ℕ) (g : F.absGal)
    (y : ↥(cyclotome (fixedFieldLocalField F S hopen).absGal m)) :
    powTorsionCongr (topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)) m
        (cyclotomeEquiv (absGalConjCME F S hopen g) m y)
      = powTorsionCongr (abelianGalConj F S hopen g) m
          (powTorsionCongr (topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen))
            m y) := by
  refine Subtype.ext ?_
  rw [powTorsionCongr_coe, powTorsionCongr_coe, powTorsionCongr_coe, cyclotomeEquiv_coe]
  show _ = topAbelianizationEquivAbelianGal _
      ((topAbelianizationCME (absGalConjCME F S hopen g)).toMulEquiv
        ((topAbelianizationEquivAbelianGal _).symm
          (topAbelianizationEquivAbelianGal _ (y : TopologicalAbelianization _))))
  rw [MulEquiv.symm_apply_apply]
  rfl

/-- `L_S ↪ K̄`(環準同型として)。★Λ9 の `torsionEquivRootsOfUnity_cyclotomic` が
要求する `ι` である。 -/
noncomputable def fixedFieldInclusion (F : PAdicLocalField p) (S : Subgroup F.absGal)
    (hopen : IsOpen (S : Set F.absGal)) :
    (fixedFieldLocalField F S hopen).carrier →+* F.closure :=
  ((IntermediateField.fixedField S).val : _ →ₐ[F.carrier] F.closure).toRingHom

/-- ★`ι ∘ σ_g = g ∘ ι` —— `σ_g` は `g` の制限なのだから定義そのもの。 -/
theorem fixedFieldInclusion_comm (F : PAdicLocalField p) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (y : (fixedFieldLocalField F S hopen).carrier) :
    fixedFieldInclusion F S hopen ((fixedFieldAut S g).toRingEquiv y)
      = g (fixedFieldInclusion F S hopen y) := rfl

/-! ## 6. ★★★★★残った壁 —— Artin 写像の単数部分の `σ`-同変性 -/

def ArtinUnitEquivariance.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★★★残った壁**(Λ12 の 5 形の、Λ9 の座標での言い換え)。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

「`μ_{p^n} ⊆ L := L_S` なる開正規 `S ⊴ Γ_F` について、局所類体論の分解
`Gal(L^{ab}/L) ≅ 𝒪_L^× × Ẑ` を、`τ ↦ g τ g^{-1}` の `𝒪_L^×` 成分が `σ_g` になるように
取れる。」

★★**`Ẑ` 成分(Frobenius 部分)は条件に現れない。** 壁は `𝒪_L^×` 成分だけである
——これは Λ12 の申し送りを本ファイルで**測って**分かったことである。

★★これは局所類体論の「Artin 写像の関手性」
`Art_L(σ_g u) = g · Art_L(u) · g^{-1}` の、`p^n` 捩れの上での姿そのものである。 -/
def ArtinUnitEquivariance (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal)
    (hopen : IsOpen (S : Set F.absGal)), S ≤ muFixer F (p ^ n) → ∀ g : F.absGal,
    ∃ E : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
        ≃* ((𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ × ZHat),
      ∀ x, x ^ (p ^ n) = 1 →
        ((unitsToField (fixedFieldLocalField F S hopen)
              (E (@abelianGalConj p _ F S hS hopen g x)).1 :
            ((fixedFieldLocalField F S hopen).carrier)ˣ)
          : (fixedFieldLocalField F S hopen).carrier)
          = fixedFieldAut S g
              ((unitsToField (fixedFieldLocalField F S hopen) (E x).1 :
                  ((fixedFieldLocalField F S hopen).carrier)ˣ)
                : (fixedFieldLocalField F S hopen).carrier)

/-! ## 7. 壁 ⇒ Λ11 の壁(★Λ9 の `torsionEquivRootsOfUnity_cyclotomic` に代入するだけ) -/

/-- ★`Ψ_g` は `tors_{p^n}(Gal(L^{ab}/L))` の上で `χ_n(g)` 乗になる。

★Λ9 の `torsionEquivRootsOfUnity_cyclotomic` に、壁の仮定 `hE` をそのまま代入した。 -/
theorem powTorsionCongr_abelianGalConj_eq_pow
    (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (E : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
        ≃* ((𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ × ZHat))
    (hE : ∀ x, x ^ (p ^ n) = 1 →
        ((unitsToField (fixedFieldLocalField F S hopen)
              (E (abelianGalConj F S hopen g x)).1 :
            ((fixedFieldLocalField F S hopen).carrier)ˣ)
          : (fixedFieldLocalField F S hopen).carrier)
          = fixedFieldAut S g
              ((unitsToField (fixedFieldLocalField F S hopen) (E x).1 :
                  ((fixedFieldLocalField F S hopen).carrier)ˣ)
                : (fixedFieldLocalField F S hopen).carrier))
    (z : ↥(powTorsion
      Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
      (p ^ n))) :
    powTorsionCongr (abelianGalConj F S hopen g) (p ^ n) z
      = z ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  refine (torsionEquivRootsOfUnity (fixedFieldLocalField F S hopen) E
    (pn_ne_zero p n)).injective ?_
  rw [torsionEquivRootsOfUnity_cyclotomic (fixedFieldLocalField F S hopen) F
      (fixedFieldInclusion F S hopen) g (fixedFieldAut S g).toRingEquiv
      (fixedFieldInclusion_comm F S hopen g) E (abelianGalConj F S hopen g) hE, map_pow]

/-- ★`Γ_{L_S}` の側の共役も `χ_n(g)` 乗になる(§1 の抽象核 A1 で輸送)。 -/
theorem cyclotomeEquiv_absGalConjCME_eq_pow
    (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (E : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
        ≃* ((𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ × ZHat))
    (hE : ∀ x, x ^ (p ^ n) = 1 →
        ((unitsToField (fixedFieldLocalField F S hopen)
              (E (abelianGalConj F S hopen g x)).1 :
            ((fixedFieldLocalField F S hopen).carrier)ˣ)
          : (fixedFieldLocalField F S hopen).carrier)
          = fixedFieldAut S g
              ((unitsToField (fixedFieldLocalField F S hopen) (E x).1 :
                  ((fixedFieldLocalField F S hopen).carrier)ˣ)
                : (fixedFieldLocalField F S hopen).carrier))
    (y : ↥(cyclotome (fixedFieldLocalField F S hopen).absGal (p ^ n))) :
    cyclotomeEquiv (absGalConjCME F S hopen g) (p ^ n) y
      = y ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) :=
  eq_pow_of_transport
    (powTorsionCongr (topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)) (p ^ n))
    (cyclotomeEquiv (absGalConjCME F S hopen g) (p ^ n))
    (powTorsionCongr (abelianGalConj F S hopen g) (p ^ n)) _
    (powTorsionCongr_abelianGalConj F S hopen (p ^ n) g)
    (powTorsionCongr_abelianGalConj_eq_pow F n S hopen g E hE) y

def cyclotomeConjIsCyclotomic_of_artinUnitEquivariance.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★壁 ⇒ Λ11 の壁**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

★Λ12 の `cyclotomeEquiv_absGalConjCME`(`Γ_{L_S}` と `↥S` の同一視)を挟むだけ。 -/
theorem cyclotomeConjIsCyclotomic_of_artinUnitEquivariance
    (h : ArtinUnitEquivariance p) : CyclotomeConjIsCyclotomic p := by
  intro F n S hS hopen hle g x
  obtain ⟨E, hE⟩ := h F n S hS hopen hle g
  obtain ⟨y, rfl⟩ := (@cyclotomeEquiv _ _ _ _ _ _ _ _ (absGalFixedFieldCME F S hopen)
    (p ^ n)).surjective x
  rw [← @cyclotomeEquiv_absGalConjCME p _ F (p ^ n) S hS hopen g y,
    @cyclotomeEquiv_absGalConjCME_eq_pow p _ F n S hS hopen g E hE y, map_pow]

/-! ## 8. ★★逆向き —— 壁は Λ11 の壁と**同値**である

★これで「言い換えで主張を強くも弱くもしていない」ことが確かめられる。 -/

/-- ★Λ11 の壁 ⇒ `Γ_{L_S}` の側の共役が `χ_n(g)` 乗。 -/
theorem cyclotomeEquiv_absGalConjCME_eq_pow_of_wall (h : CyclotomeConjIsCyclotomic p)
    (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal)
    (hopen : IsOpen (S : Set F.absGal)) (hle : S ≤ muFixer F (p ^ n)) (g : F.absGal)
    (y : ↥(cyclotome (fixedFieldLocalField F S hopen).absGal (p ^ n))) :
    cyclotomeEquiv (@absGalConjCME p _ F S hS hopen g) (p ^ n) y
      = y ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  haveI := hS
  refine (@cyclotomeEquiv _ _ _ _ _ _ _ _ (absGalFixedFieldCME F S hopen)
    (p ^ n)).injective ?_
  rw [@cyclotomeEquiv_absGalConjCME p _ F (p ^ n) S hS hopen g y,
    h F n S hS hopen hle g _, map_pow]

/-- ★★Λ11 の壁 ⇒ `Ψ_g` は捩れの上で `χ_n(g)` 乗(**元の言葉で**)。

★★部分型を経由せずに `Gal(L^{ab}/L)` の元の等式として書いた形。
これが §8 の残りで実際に使う姿である。 -/
theorem abelianGalConj_eq_pow_of_wall (h : CyclotomeConjIsCyclotomic p)
    (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal)
    (hopen : IsOpen (S : Set F.absGal)) (hle : S ≤ muFixer F (p ^ n)) (g : F.absGal)
    (x : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier))
    (hx : x ^ (p ^ n) = 1) :
    @abelianGalConj p _ F S hS hopen g x
      = x ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  haveI := hS
  have hw : ((topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)).symm x)
      ^ (p ^ n) = 1 := by
    refine (topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)).injective ?_
    rw [map_pow, MulEquiv.apply_symm_apply, hx, map_one]
  have hy := cyclotomeEquiv_absGalConjCME_eq_pow_of_wall h F n S hS hopen hle g
    ⟨(topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)).symm x, hw⟩
  have hy' : topAbelianizationCME (@absGalConjCME p _ F S hS hopen g)
        ((topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)).symm x)
      = ((topAbelianizationEquivAbelianGal (fixedFieldLocalField F S hopen)).symm x)
          ^ ((PadicInt.toZModPow n
            ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) :=
    congrArg Subtype.val hy
  rw [abelianGalConj_apply, hy', map_pow, MulEquiv.apply_symm_apply]

/-- ★★★**逆向き —— Λ11 の壁 ⇒ 本ファイルの壁**。

`μ_{p^n}(L_S)` の上では `σ_g` は `χ_n(g)` 倍である
(Λ12 の `rootsOfUnityFixedFieldConj_eq_pow`、★ここで `S ≤ muFixer` を使う)。
左辺も `Ψ_g` が `χ_n(g)` 乗であること(`abelianGalConj_eq_pow_of_wall`)から
`χ_n(g)` 乗になる。★**両辺が同じ指数になるので一致する。** -/
theorem artinUnitEquivariance_of_cyclotomeConjIsCyclotomic
    (h : CyclotomeConjIsCyclotomic p) : ArtinUnitEquivariance p := by
  intro F n S hS hopen hle g
  haveI := hS
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  obtain ⟨Φ⟩ := nonempty_abelianGalContinuousEquivUnitsZHat (fixedFieldLocalField F S hopen)
  refine ⟨Φ.toMulEquiv, fun x hx => ?_⟩
  have hfst : ∀ (a : (𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ × ZHat) (k : ℕ),
      (a ^ k).1 = a.1 ^ k := fun _ _ => rfl
  have hΨ : @abelianGalConj p _ F S hS hopen g x
      = x ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) :=
    abelianGalConj_eq_pow_of_wall h F n S hS hopen hle g x hx
  have hΦpow : ∀ k : ℕ, Φ.toMulEquiv (x ^ k) = (Φ.toMulEquiv x) ^ k := fun k => map_pow _ _ _
  have hleft : (unitsToField (fixedFieldLocalField F S hopen)
        (Φ.toMulEquiv (@abelianGalConj p _ F S hS hopen g x)).1
      : ((fixedFieldLocalField F S hopen).carrier)ˣ)
      = (unitsToField (fixedFieldLocalField F S hopen) (Φ.toMulEquiv x).1)
          ^ ((PadicInt.toZModPow n
            ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
    rw [hΨ, hΦpow, hfst, map_pow]
  have hu : (unitsToField (fixedFieldLocalField F S hopen) (Φ.toMulEquiv x).1) ^ (p ^ n) = 1 := by
    have h1 : (Φ.toMulEquiv x) ^ (p ^ n) = 1 := by rw [← hΦpow, hx, map_one]
    have h2 : ((Φ.toMulEquiv x).1) ^ (p ^ n) = 1 := by rw [← hfst, h1]; rfl
    rw [← map_pow, h2, map_one]
  have humem : (unitsToField (fixedFieldLocalField F S hopen) (Φ.toMulEquiv x).1)
      ∈ rootsOfUnity (p ^ n) ↥(IntermediateField.fixedField S) := (mem_rootsOfUnity _ _).mpr hu
  have hright := rootsOfUnityFixedFieldConj_eq_pow F n (S := S) hle g
    ⟨unitsToField (fixedFieldLocalField F S hopen) (Φ.toMulEquiv x).1, humem⟩
  have e1 : (((rootsOfUnityFixedFieldConj F (p ^ n) S g
        ⟨unitsToField (fixedFieldLocalField F S hopen) (Φ.toMulEquiv x).1, humem⟩ :
        ↥(rootsOfUnity (p ^ n) ↥(IntermediateField.fixedField S)))
        : (↥(IntermediateField.fixedField S))ˣ) : ↥(IntermediateField.fixedField S))
      = fixedFieldAut S g
          ((unitsToField (fixedFieldLocalField F S hopen) (Φ.toMulEquiv x).1
            : ((fixedFieldLocalField F S hopen).carrier)ˣ)
            : (fixedFieldLocalField F S hopen).carrier) :=
    rootsOfUnityFixedFieldConj_coe F (p ^ n) S g _
  rw [hleft, ← e1, hright, SubmonoidClass.coe_pow, Units.val_pow_eq_pow_val,
    Units.val_pow_eq_pow_val]
  rfl

/-- **★★★★★★★★本ファイルの壁は Λ11 の壁と同値**。

★言い換えで主張を強くも弱くもしていない。★これで Λ12 の 5 形と合わせて
**6 つの同値な言い方**が揃った。 -/
theorem artinUnitEquivariance_iff_cyclotomeConj :
    ArtinUnitEquivariance p ↔ CyclotomeConjIsCyclotomic p :=
  ⟨cyclotomeConjIsCyclotomic_of_artinUnitEquivariance,
    artinUnitEquivariance_of_cyclotomeConjIsCyclotomic⟩

/-! ## 9. Λ12 の 5 形・Proposition 1.1 への配線 -/

/-- ★本ファイルの壁 ⇒ Λ12 の `ArtinEquivarianceLocalField`。 -/
theorem artinEquivarianceLocalField_of_artinUnitEquivariance
    (h : ArtinUnitEquivariance p) : ArtinEquivarianceLocalField p :=
  artinEquivarianceLocalField_of_fixedField
    (artinEquivarianceFixedField_iff_cyclotomeConj.mpr
      (cyclotomeConjIsCyclotomic_of_artinUnitEquivariance h))

/-- ★本ファイルの壁 ⇒ Λ12 の `ArtinEquivariance`。 -/
theorem artinEquivariance_of_artinUnitEquivariance
    (h : ArtinUnitEquivariance p) : ArtinEquivariance p :=
  artinEquivariance_iff_cyclotomeConj.mpr
    (cyclotomeConjIsCyclotomic_of_artinUnitEquivariance h)

def cyclotomicCharacter_recoverable_of_artinUnitEquivariance.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★★★★★★★★★★★★[pGC] Proposition 1.1 の現在地**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

★★**Proposition 1.1 は「Artin 写像の `𝒪_L^×` 成分が `σ_g` と同変」1 点だけに依存している。**
★`Skeleton/PGC/Section1.lean` の `sorry` は**まだ埋まっていない**。 -/
theorem cyclotomicCharacter_recoverable_of_artinUnitEquivariance
    (h : ArtinUnitEquivariance p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal :=
  cyclotomicCharacter_recoverable_of_cyclotomeConj
    (cyclotomeConjIsCyclotomic_of_artinUnitEquivariance h)

/-! ## 10. ★退化の自己検査 —— 壁の易しい半分(`g ∈ S`)は実際に成り立つ -/

/-- ★**`g ∈ S` なら `Ψ_g` は恒等**。

`g ∈ S` による共役は `Γ_{L_S}^{ab}` の上で恒等(Λ2 の `topAbelianizationCME_conj_self`)。 -/
theorem abelianGalConj_of_mem (F : PAdicLocalField p) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) {g : F.absGal} (hg : g ∈ S)
    (x : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)) :
    abelianGalConj F S hopen g x = x := by
  rw [abelianGalConj_apply]
  have hstep : ∀ w : TopologicalAbelianization (fixedFieldLocalField F S hopen).absGal,
      topAbelianizationCME (absGalConjCME F S hopen g) w = w := by
    intro w
    have h1 : topAbelianizationCME (absGalConjCME F S hopen g) w
        = topAbelianizationCME (absGalFixedFieldCME F S hopen).symm
            (topAbelianizationCME
              (conjSubgroupCME (H := S) ((⟨g, hg⟩ : ↥S) : F.absGal))
              (topAbelianizationCME (absGalFixedFieldCME F S hopen) w)) := by
      rw [topAbelianizationCME_trans, topAbelianizationCME_trans]
      rfl
    rw [h1, topAbelianizationCME_conj_self ⟨g, hg⟩, topAbelianizationCME_trans]
    exact (topAbelianizationCME_congr
      (α := (absGalFixedFieldCME F S hopen).trans (absGalFixedFieldCME F S hopen).symm)
      (β := ContinuousMulEquiv.refl _)
      (fun a => (absGalFixedFieldCME F S hopen).symm_apply_apply a) w).trans
      (topAbelianizationCME_refl w)
  rw [hstep, MulEquiv.apply_symm_apply]

/-- ★★**壁の易しい半分** —— `g ∈ S` なら `ArtinUnitEquivariance` の条件は
**どの `E` についても**成り立つ。

★左辺は `abelianGalConj_of_mem`、右辺は `fixedFieldAut_of_mem` で、
**独立の理由で**両辺が恒等になる。★したがって壁の内容は
`Γ_F ⧸ S = Gal(L_S/F)` の上にある。 -/
theorem artinUnitEquivariance_condition_of_mem (F : PAdicLocalField p) (n : ℕ)
    (S : Subgroup F.absGal) [S.Normal] (hopen : IsOpen (S : Set F.absGal))
    {g : F.absGal} (hg : g ∈ S)
    (E : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
        ≃* ((𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ × ZHat))
    (x : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier))
    (_hx : x ^ (p ^ n) = 1) :
    ((unitsToField (fixedFieldLocalField F S hopen)
          (E (abelianGalConj F S hopen g x)).1 :
        ((fixedFieldLocalField F S hopen).carrier)ˣ)
      : (fixedFieldLocalField F S hopen).carrier)
      = fixedFieldAut S g
          ((unitsToField (fixedFieldLocalField F S hopen) (E x).1 :
              ((fixedFieldLocalField F S hopen).carrier)ˣ)
            : (fixedFieldLocalField F S hopen).carrier) := by
  rw [abelianGalConj_of_mem F S hopen hg x, fixedFieldAut_of_mem S hg]

/-- ★**易しい半分の存在形** —— `g ∈ S` については壁が実際に埋まっている。 -/
theorem exists_artinUnitEquivariance_of_mem (F : PAdicLocalField p) (n : ℕ)
    (S : Subgroup F.absGal) [S.Normal] (hopen : IsOpen (S : Set F.absGal))
    {g : F.absGal} (hg : g ∈ S) :
    ∃ E : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
        ≃* ((𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ × ZHat),
      ∀ x, x ^ (p ^ n) = 1 →
        ((unitsToField (fixedFieldLocalField F S hopen)
              (E (abelianGalConj F S hopen g x)).1 :
            ((fixedFieldLocalField F S hopen).carrier)ˣ)
          : (fixedFieldLocalField F S hopen).carrier)
          = fixedFieldAut S g
              ((unitsToField (fixedFieldLocalField F S hopen) (E x).1 :
                  ((fixedFieldLocalField F S hopen).carrier)ˣ)
                : (fixedFieldLocalField F S hopen).carrier) := by
  obtain ⟨Φ⟩ := nonempty_abelianGalContinuousEquivUnitsZHat (fixedFieldLocalField F S hopen)
  exact ⟨Φ.toMulEquiv, fun x hx =>
    artinUnitEquivariance_condition_of_mem F n S hopen hg Φ.toMulEquiv x hx⟩

/-! ## 11. ★抽象核 D —— Weierstrass 標準分解の自然性

★★「点で評価する層」を木の捩れ点(＝`iteratedLubinTateDistinguished` の根)に
届かせるための最後の一般論。★この節にも体・Galois・分岐・Lubin-Tate の語彙は
1 つも出てこない —— 一般の可換局所環である。

★mathlib の `PowerSeries.IsWeierstrassFactorization.unique`(分解の一意性)に
乗せるだけで出る。 -/

section WeierstrassTwist

variable {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [IsLocalRing B]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]

omit [IsAdicComplete (IsLocalRing.maximalIdeal A) A] in
/-- ★抽象核 D1 —— 局所環の同型は極大イデアルを極大イデアルに写す(元ごとの形)。 -/
theorem map_mem_maximalIdeal (φ : A ≃+* B) {a : A} (ha : a ∈ IsLocalRing.maximalIdeal A) :
    φ a ∈ IsLocalRing.maximalIdeal B := by
  rw [IsLocalRing.mem_maximalIdeal] at ha ⊢
  intro hu
  exact ha (by simpa using hu.map φ.symm)

omit [IsLocalRing A] [IsLocalRing B] [IsAdicComplete (IsLocalRing.maximalIdeal A) A] in
/-- ★多項式から冪級数への包含は係数写像と可換。 -/
theorem coe_polynomial_map (φ : A →+* B) (P : Polynomial A) :
    ((P.map φ : Polynomial B) : PowerSeries B) = PowerSeries.map φ (P : PowerSeries A) := by
  ext n
  rw [Polynomial.coeff_coe, Polynomial.coeff_map, PowerSeries.coeff_map, Polynomial.coeff_coe]

omit [IsAdicComplete (IsLocalRing.maximalIdeal A) A] in
/-- ★抽象核 D2 —— distinguished 多項式は係数をひねっても distinguished。 -/
theorem isDistinguishedAt_map (φ : A ≃+* B) {P : Polynomial A}
    (hP : P.IsDistinguishedAt (IsLocalRing.maximalIdeal A)) :
    (P.map (φ : A →+* B)).IsDistinguishedAt (IsLocalRing.maximalIdeal B) := by
  refine ⟨⟨fun {n} hn => ?_⟩, hP.monic.map _⟩
  rw [Polynomial.coeff_map]
  refine map_mem_maximalIdeal φ (hP.mem ?_)
  rwa [hP.monic.natDegree_map] at hn

/-- **★★抽象核 D3 —— Weierstrass 標準分解は係数のひねりと可換**:

  `(g^φ)` の distinguished 部分 `=` (`g` の distinguished 部分)`^φ`。

★分解の**一意性**(mathlib `IsWeierstrassFactorization.unique`)に
「`g^φ = (D^φ)·(U^φ)` は Weierstrass 分解である」を代入するだけ。 -/
theorem map_weierstrassDistinguished [IsAdicComplete (IsLocalRing.maximalIdeal B) B]
    (φ : A ≃+* B) {g : PowerSeries A}
    (hg : PowerSeries.map (IsLocalRing.residue A) g ≠ 0)
    (hg' : PowerSeries.map (IsLocalRing.residue B) (PowerSeries.map (φ : A →+* B) g) ≠ 0) :
    (PowerSeries.map (φ : A →+* B) g).weierstrassDistinguished hg'
      = Polynomial.map (φ : A →+* B) (g.weierstrassDistinguished hg) := by
  have H : (PowerSeries.map (φ : A →+* B) g).IsWeierstrassFactorization
      (Polynomial.map (φ : A →+* B) (g.weierstrassDistinguished hg))
      (PowerSeries.map (φ : A →+* B) (g.weierstrassUnit hg)) := by
    refine ⟨isDistinguishedAt_map φ (PowerSeries.isDistinguishedAt_weierstrassDistinguished hg),
      (PowerSeries.isUnit_weierstrassUnit hg).map (PowerSeries.map (φ : A →+* B)), ?_⟩
    rw [coe_polynomial_map, ← map_mul]
    exact congrArg (PowerSeries.map (φ : A →+* B))
      (PowerSeries.eq_weierstrassDistinguished_mul_weierstrassUnit hg)
  exact (H.unique hg').1.symm

end WeierstrassTwist

section WeierstrassCongr

variable {A : Type*} [CommRing A] [IsLocalRing A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]

/-- ★`weierstrassDistinguished` は第 1 引数にしか依らない(証明部分は無関係)。 -/
theorem weierstrassDistinguished_congr {g g' : PowerSeries A} (h : g = g')
    (hg : PowerSeries.map (IsLocalRing.residue A) g ≠ 0)
    (hg' : PowerSeries.map (IsLocalRing.residue A) g' ≠ 0) :
    g.weierstrassDistinguished hg = g'.weierstrassDistinguished hg' := by
  subst h; rfl

end WeierstrassCongr

/-! ## 12. 具体層 3 —— `D_n^{f^φ} = (D_n^f)^φ`

★★木の捩れ点は `iteratedLubinTateDistinguished`(`[π^n]_f` の Weierstrass
distinguished 部分)の根として定義されているので、「`σ` が捩れ点を運ぶ」ためには
**この自然性**が要る。★在庫 `map_iteratedLubinTate`(`DworkThetaEval.lean`、
係数写像は自己合成と可換)と §11 の抽象核 D3 を継ぐだけ。 -/

section IteratedTwist

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp]
    [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)

/-- **★★★`D_n` の自然性** —— `f` の係数を `φ` でひねると `D_n` もひねられる。 -/
theorem map_iteratedLubinTateDistinguished (φ : A ≃+* A)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    iteratedLubinTateDistinguished hq (maximalIdeal_eq_span_map φ hπmax)
        (by simpa using fun h => hπne0 ((φ : A ≃+* A).injective (by rw [h, map_zero])))
        (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
        (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n
      = Polynomial.map (φ : A →+* A)
          (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by simpa using hf0
  have hnat := map_iteratedLubinTate (φ : A →+* A) hf0c n
  have h1 : iteratedLubinTateDistinguished hq (maximalIdeal_eq_span_map φ hπmax)
        (by simpa using fun h => hπne0 ((φ : A ≃+* A).injective (by rw [h, map_zero])))
        (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
        (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n
      = (PowerSeries.map (φ : A →+* A) (iteratedLubinTate f n)).weierstrassDistinguished
          (by rw [hnat]
              exact iteratedLubinTate_map_residue_ne_zero hq (maximalIdeal_eq_span_map φ hπmax)
                (by simpa using fun h => hπne0 ((φ : A ≃+* A).injective (by rw [h, map_zero])))
                (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
                (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n) :=
    weierstrassDistinguished_congr hnat.symm _ _
  rw [h1, map_weierstrassDistinguished φ
    (iteratedLubinTate_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n)]
  rfl

end IteratedTwist

/-! ## 13. ★★★★「点で評価する層」の到達点 —— `σ` は `Λ_{f,n}` を `Λ_{f^σ,n}` に運ぶ

★★これが Λ12 が名指しした入力 (b) である。 -/

open scoped Classical in
/-- **★★★★★`σ` は `f`-捩れ点を `f^σ`-捩れ点に写す**。

`Φ : K̄ → K̄` が `φ : 𝒪_K ≃ 𝒪_K` について半線型(`Φ ∘ ι = ι ∘ φ`)なら、
`Λ_{f,n}` の元は `Λ_{f^φ,n}` に写る。

★★これが Λ12 の申し送り (b)「`σ` が捩れ点 `μ_{f,m} → μ_{f^σ,m}` を運ぶ
**点で評価する層**」の中身である。3 段の合成:
§3 の抽象核 C2(根は根に写る)+ §12 の `D_n` の自然性 +
`iteratedLubinTateTorsionPoints` の定義(`D_n` の根の `Finset`)。

★**`Φ` に `K`-線型性を仮定していない**(仮定すると `φ = id` になって主張が空になる)。
実際に代入するのは `Φ := g ∈ Γ_F`、`φ := fixedFieldIntegerAut S g` である。 -/
theorem map_mem_iteratedLubinTateTorsionPoints (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (φ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier])
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (Φ : K.closure →+* K.closure)
    (hcompat : ∀ a : 𝒪[K.carrier], Φ (algebraMap 𝒪[K.carrier] K.closure a)
      = algebraMap 𝒪[K.carrier] K.closure (φ a))
    {x : K.closure}
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    Φ x ∈ iteratedLubinTateTorsionPoints K hq (maximalIdeal_eq_span_map φ hπmax)
      (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
      (PowerSeries.map (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f) (coeff_zero_map_twist _ hf0)
      (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n := by
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hx ⊢
  obtain ⟨-, hroot⟩ := hx
  refine ⟨?_, ?_⟩
  · refine Polynomial.Monic.ne_zero ?_
    exact ((isDistinguishedAt_iteratedLubinTateDistinguished hq (maximalIdeal_eq_span_map φ hπmax)
      (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
      (PowerSeries.map (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f) (coeff_zero_map_twist _ hf0)
      (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n).monic).map _
  · rw [Polynomial.IsRoot.def, Polynomial.eval_map] at hroot ⊢
    rw [map_iteratedLubinTateDistinguished hq φ hπmax hπne0 f hf0 hf1 hf n]
    exact eval₂_map_twist_eq_zero (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) Φ
      (algebraMap 𝒪[K.carrier] K.closure) hcompat _ hroot

end ABC3.Found.PGC
