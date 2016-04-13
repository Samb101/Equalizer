%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% AUTO GENERE PAR MATLAB NON MODIFIE %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = Projet(varargin)

    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Projet_OpeningFcn, ...
                       'gui_OutputFcn',  @Projet_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end



function Projet_OpeningFcn(hObject, eventdata, handles, varargin)

    handles.output = hObject;
    
    guidata(hObject, handles);
    
function varargout = Projet_OutputFcn(hObject, eventdata, handles) 

    varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nos FONCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Ouverture, lecture, tracage des courbes %%%%%%%%%%%%%%%%%%%

function openFilebutton_Callback(hObject, eventdata, handles)
    % fonction d'ouverture et de mise en route de la lecture et des
    % graphiques
    handles = guidata(hObject);
   
    %ouverture du fichier
    [handles.file_name] = uigetfile('*.m4a','File Selector');
    [handles.audio, handles.fr] = audioread(handles.file_name);
    handles.origin_tf = fft(handles.audio);

    % declaration de l'audioplayer et des handles de 'sauvegarde'
    handles.playlist = audioplayer(handles.audio, handles.fr);
    handles.current_audio = handles.audio; 
    handles.current_tf = handles.origin_tf;
    
    % affiche le titre de la chanson et sa frequence
    set(handles.uipanel5, 'Visible', 'on');
    set(handles.path, 'String', handles.file_name, 'FontWeight', 'Bold');
    set(handles.frequency, 'String', strcat(num2str(handles.fr), 'Hz'));
    
    play(handles.playlist);

    set(handles.playbutton, 'String', 'Pause');
    
    guidata(hObject, handles);
    
    clearaxes(handles);
    sync(handles);

    
function clearaxes(handles)
    % courte fonction pour nettoyer les graphiques
    % en cas de modification de fr?quence on peut revenir 
    % sur une plage deja jouee donc il faut clean les axes
    cla(handles.temp,'reset');
    cla(handles.freq,'reset');
        
function sync(handles)
    % Fonction qui trace les courbes
    
    % declaration des lignes anim?es
    axes(handles.temp);
    h_original_temp = animatedline('Color','b');
    h_current_temp = animatedline('Color','r');
    
    axes(handles.freq);
    h_original_freq = animatedline;
    h_current_freq = animatedline('Color','g');
 
    % boucle de tracage des courbes
    for k=1:length(handles.playlist.TotalSamples)
        while isplaying(handles.playlist) == true  
                
            % affiche le temps en lettres et sur le slider
            handles.audiotime.String = strcat('Time :', num2str(ceil(handles.playlist.CurrentSample/handles.playlist.SampleRate)), 's');
            handles.timeslider.Value = handles.playlist.CurrentSample/handles.playlist.TotalSamples;
                
            
            % incremente les axes tous les multiples de la frequence * 2,
            % soit 2 secondes
            if  mod(handles.fr*2,k) == 0
                axis(handles.temp, [((handles.playlist.CurrentSample-handles.fr*2)/(handles.playlist.SampleRate')) ((handles.playlist.CurrentSample'+handles.fr*2)/(handles.playlist.SampleRate')) -1 1]);
                axis(handles.freq, [((get(handles.playlist, 'CurrentSample')-handles.fr*2)/(get(handles.playlist,'SampleRate'))) ((get(handles.playlist, 'CurrentSample')+handles.fr*2)/(get(handles.playlist,'SampleRate'))) 0 1000]);
            end
            
            addpoints(h_current_temp, handles.playlist.CurrentSample/handles.playlist.SampleRate, handles.current_audio(get(handles.playlist, 'CurrentSample'))');
            if handles.original_checkbox.Value == 1 % si on veut comparer
                addpoints(h_original_temp, handles.playlist.CurrentSample/handles.fr, handles.audio(get(handles.playlist, 'CurrentSample'))');
            end
            addpoints(h_current_freq, handles.playlist.CurrentSample/handles.playlist.SampleRate, fftshift(abs(handles.current_tf(get(handles.playlist, 'CurrentSample'))')));
            if handles.original_freq_checkbox.Value == 1
                addpoints(h_original_freq, handles.playlist.CurrentSample/handles.fr, abs(handles.origin_tf(get(handles.playlist, 'CurrentSample'))'));
            end
            drawnow;

        end
    end

  
function playbutton_Callback(hObject, eventdata, handles)
    % bouton play/pause/resume tout en un
    if isplaying(handles.playlist) && strcmp ('Pause', get(handles.playbutton, 'String'))

        pause(handles.playlist);
        set(handles.playbutton,'String','Resume');
        return;

    else
      
        resume(handles.playlist);
        set(handles.playbutton,'String','Pause');
        sync(handles);
        return;
      
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Fonctions temporelles  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%    Utilitaires  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timeslider_Callback(hObject, eventdata, handles)
    % Permet d'avancer/reculer dans la chanson
    % marche mal car la valeur est en m?me temps set par la fonction sync()
    % il faut spammer le clic mais cela fait parfois plante le logiciel
    if isplaying(handles.playlist)
       
        pause(handles.playlist);
        play(handles.playlist, ceil(get(handles.timeslider,'Value')*get(handles.playlist,'TotalSamples')));

        guidata(hObject,handles);

        clearaxes(handles);
        sync(handles);
    end


 function freqslider_Callback(hObject, eventdata, handles)
    % modifie la frequence de l'audio courant
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        a = handles.playlist.CurrentSample;
        b = handles.playlist.SampleRate;
        pause(handles.playlist);
        
        handles.playlist.SampleRate = handles.fr*get(handles.freqslider, 'Value')*2;
        set(handles.frequency, 'String', b);
        
        guidata(hObject, handles);
        
        play(handles.playlist, a);
        
        clearaxes(handles);
        sync(handles);
    end

function volumeslider_Callback(hObject, eventdata, handles)
    % modifie le volume de l'audioplayer entre x0 et x2
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        a = handles.playlist.CurrentSample;
        b = handles.playlist.SampleRate;
        pause(handles.playlist);

        y = handles.current_audio * handles.volumeslider.Value * 2;    
        
        % mise a jour
        handles.current_audio = y;
        handles.current_tf = fft(handles.current_audio);
        handles.playlist = audioplayer(handles.current_audio,b);
       
        guidata(hObject,handles); 
        
        play(handles.playlist, a);
        
        sync(handles);
        
    end

function backtooriginbutton_Callback(hObject, eventdata, handles)
    % permet de revenir au son d'origine
    if isplaying(handles.playlist) == true 
       pause(handles.playlist); 
    end
    handles.current_audio = handles.audio;
    handles.current_tf = fft(handles.current_audio);
    handles.playlist = audioplayer(handles.current_audio, handles.fr);
    guidata(hObject,handles);
    clearaxes(handles);
    play(handles.playlist);
    sync(handles);
    

%%%%%%%%%%%%%%%%%%%%%%    Effets temporels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function revertbutton_Callback(hObject, eventdata, handles)
    % lis la matrice inversee
    if isplaying(handles.playlist)
        
        a = handles.playlist.CurrentSample;
        b = handles.playlist.SampleRate;
        pause(handles.playlist);
        
        handles.playlist = audioplayer(flipud(handles.current_audio), b);
        % flipud inverse la matrice
        guidata(hObject,handles);

        play(handles.playlist, a);

        clearaxes(handles);
        sync(handles);
    end

    

function echoslider_Callback(hObject, eventdata, handles)
    % ajoute un effet d'echo a l'audio courante                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
     if isplaying(handles.playlist)
        % les donnees de lecture courante
         a = handles.playlist.CurrentSample; 
         b = handles.playlist.SampleRate;
         pause(handles.playlist);

        
         x = handles.current_audio;
         N = handles.fr*ceil((handles.echoslider.Value*5));
         % l'ecart d'echo est defini par la valeur du slider
         
         for n = N+1 : length(handles.current_audio)
              x(n,1) = handles.current_audio(n,1) + handles.current_audio(n-N,1);
              x(n,2) = handles.current_audio(n,2) + handles.current_audio(n-N,2);
         end
         
         % mise a jour
         handles.current_audio = x;
         handles.current_tf = fft(handles.current_audio);

         handles.playlist= audioplayer(handles.current_audio, b);
         play(handles.playlist,a);

                  
         guidata(hObject,handles);

         clearaxes(handles);
         sync(handles);
     end
        
function tremolobutton_Callback(hObject, eventdata, handles)
    % effet tremolo sur l'audio courant
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        a = handles.playlist.CurrentSample;
        pause(handles.playlist);
        x = handles.current_audio;
        Fe = handles.playlist.SampleRate;

        Fm=5;
        A=0.5;

        i1 = 1:length(x(:,1));
        i2 = 1:length(x(:,2));
        h1 = (1 + A*sin(2*pi*i1*(Fm/Fe)))';
        h2 = (1 + A*sin(2*pi*i2*(Fm/Fe)))';
        y(:,1) = times(h1,x(:,1));
        y(:,2) = times(h2,x(:,2));

        
        % mise ? jour
        handles.current_audio = y;
        handles.current_tf = fft(handles.current_audio);
        handles.playlist = audioplayer(handles.current_audio,Fe);
       
        guidata(hObject,handles); 
        
        play(handles.playlist, a);
        
        sync(handles);
        
    end

            
function wahwahbutton_Callback(hObject, eventdata, handles)
    % effet wawah sur l'audio courant
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        a = handles.playlist.CurrentSample;
        pause(handles.playlist);
        x = handles.current_audio;
        Fe = handles.playlist.SampleRate;
     

        amortfact=0.05;
        minf=500;
        maxf=3000;

        Fw = 2000;

        delta = Fw/Fe;

        Fc=minf:delta:maxf;
        
        while(length(Fc) < length(x) )
            Fc= [ Fc (maxf:-delta:minf) ];
            Fc= [ Fc (minf:delta:maxf) ];
        end

        Fc = Fc(1:length(x));

        F1 = 2*sin((pi*Fc(1))/Fe);

        Q1 = 2*amortfact;
        yh=zeros(size(x)); 
        yb=zeros(size(x));
        yl=zeros(size(x));

        % le premier element est a blinder il peut etre parfois n?gatif
        yh(1) = x(1);
        yb(1) = F1*yh(1);
        yl(1) = F1*yb(1);

        for n=2:length(x),
            yh(n) = x(n) - yl(n-1) - Q1*yb(n-1);
            yb(n) = F1*yh(n) + yb(n-1);
            yl(n) = F1*yb(n) + yl(n-1);
            F1 = 2*sin((pi*Fc(n))/Fe);
        end

        %normalisation
        maxyb = max(abs(yb));
        y = yb/maxyb;

        % mise a jour
        handles.current_audio = y;
        handles.current_tf = fft(handles.current_audio);
        handles.playlist = audioplayer(handles.current_audio,Fe);
       
        guidata(hObject,handles); 
        
        play(handles.playlist,a);
        
        sync(handles);
        
    end

function vibratobutton_Callback(hObject, eventdata, handles)
    % effet vibratobutton sur l'audio courant
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        a = handles.playlist.CurrentSample;
        pause(handles.playlist);
        x = handles.current_audio;
        Fe = handles.playlist.SampleRate;

        f0=10;
        delta=0.0008;

        % on cree l'audio tampon (vecteur de zeros) et le time_delay de la
        % vibration
        y =  zeros(size(x));
        DELAY = round(delta.*Fe);
        frequence = f0/Fe; 
        long_delay=2+DELAY*3;
        DL=zeros(long_delay,1); 

        for n=1: (length(x)-1)
            MOD=sin(frequence*2*pi*n);
            zeiger=1+DELAY+DELAY*MOD;
            i=floor(zeiger);
            frac=zeiger-i;
            DL=[x(n);DL(1:long_delay-1)];
   
            y(n,1)=DL(i+1)*frac+DL(i)*(1-frac);
        end
        
        % mise a jour
        handles.current_audio = y;
        handles.current_tf = fft(handles.current_audio);
        handles.playlist = audioplayer(handles.current_audio,Fe);
        
        guidata(hObject,handles); 
        
        play(handles.playlist,a);
        
        sync(handles);

    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Fonctions frequentielles  %%%%%%%%%%%%%%%%%%%%%

function passe_bas_slider_Callback(hObject, eventdata, handles)
    % filtre passe-bas, le slider r?presente la bande passante K
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        a = handles.playlist.SampleRate;
        b = handles.playlist.CurrentSample;
        pause(handles.playlist);
        
        % analogue au TP6
        K = ceil(handles.passe_bas_slider.Value*3000000);
        Y =  zeros(size(handles.current_tf)); 
        Y(1:K, 1) = handles.current_tf(1:K, 1); 
        Y(1:K, 2) = handles.current_tf(1:K, 2); 
        
        % mise a jour
        handles.current_audio = real(ifft(Y));
        handles.current_tf = Y;
        handles.playlist = audioplayer(handles.current_audio, a);
        guidata(hObject, handles);
        
        play(handles.playlist, b);
            
        sync(handles);
        
    end
    
function highpassslider_Callback(hObject, eventdata, handles)
 % filtre passe-bas, le slider r?presente la bande passante K
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        a = handles.playlist.SampleRate;
        b = handles.playlist.CurrentSample;
        x = handles.current_audio;
        pause(handles.playlist);
        
        % analogue au TP6
        K = ceil(handles.highpassslider.Value*3000000);
        Y =  zeros(size(handles.current_tf)); 
        Y(length(x)- 1 -K : length(x) - 1, 1) = handles.current_tf(length(x) -1 -K :length(x) - 1, 1); 
        Y(length(x)- 1 -K : length(x) - 1, 2) = handles.current_tf(length(x) -1 -K :length(x) - 1, 2); 

         
        % mise a jour
        handles.current_audio = real(ifft(Y));
        handles.current_tf = Y;
        handles.playlist = audioplayer(handles.current_audio, a);
        guidata(hObject, handles);
        
        play(handles.playlist, b);
            
        sync(handles);
        
    end

function ordreslider_Callback(hObject, eventdata, handles)
    set(handles.ordrestring, 'String', strcat('Ordre', ':', num2str(ceil(handles.ordreslider.Value*5))));
    
function butterworthbutton_Callback(hObject, eventdata, handles)
   % filtre de butterworth
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        g = handles.playlist.CurrentSample;
        pause(handles.playlist);
        x = handles.current_audio;
        Fe = handles.playlist.SampleRate;
        ordre = ceil(handles.ordreslider.Value * 5); 

        % analogue au TP3
        w= 400/Fe;
        d=0.9;

        [b,a] = butter(ordre,[.15,d]);
        [c,d] = butter(ordre,w,'s');
        
        y=filter(b,a,x);

        mag=bode(c,d);

        plot(handles.bode, squeeze(mag).^2);
        xlabel(handles.bode, 'Frequence (Hz)'); 
        ylabel(handles.bode, '|H(jOmega)|^2   (Puissance de la reponse du filtre)');
        title(handles.bode, 'Diagramme de Bode');
        
        filt = tf( c , d );
        sgrid
        pzmap(handles.poleszero,filt); 
        grid on 
        axis equal
        title(handles.poleszero, 'Diagramme de poles zeros');
   
        % mise a jour
        handles.current_audio = y;
        handles.current_tf = fft(y);
        handles.playlist = audioplayer(handles.current_audio,Fe);
        
        guidata(hObject,handles); 
        
        play(handles.playlist,g);
        
        sync(handles);
    
    end
    
function tchebychevbutton_Callback(hObject, eventdata, handles)
    %filtre de tchebychev
    if isplaying(handles.playlist)
        % les donnees de lecture courante
        g = handles.playlist.CurrentSample;
        pause(handles.playlist);
        x = handles.current_audio;
        Fe = handles.playlist.SampleRate;
        ordre = ceil(handles.ordreslider.Value * 5); 
        
        % analogue au TP3
        w= 400/Fe;
        d=0.9;

        [b,a] = cheby1(ordre,2,[.15,d]);
        [c,d] = cheby1(ordre,4,w,'s');
        
        y=filter(b,a,x);

        mag=bode(c,d);        
        
        
        plot(handles.bode, squeeze(mag).^2);
        xlabel(handles.bode, 'Frequence (Hz)'); 
        ylabel(handles.bode, '|H(jOmega)|^2   (Puissance de la reponse du filtre)');
        title(handles.bode, 'Diagramme de Bode');
        
        
        filt = tf( c , d );
        sgrid
        pzmap(handles.poleszero,filt); 
        grid on 
        axis equal
        title(handles.poleszero, 'Diagramme de poles zeros');
   
        % mise a jour
        handles.current_audio = y;
        handles.current_tf = fft(y);
        handles.playlist = audioplayer(handles.current_audio,Fe);
        
        guidata(hObject,handles); 
        
        play(handles.playlist,g);
        sync(handles);

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Fonctions d'export  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function savebutton_Callback(hObject, eventdata, handles)
        % on ouvre une fenetre de l'OS pour enregistrer le fichier
    FileName = uiputfile('*.m4a', 'Save your fantastic audio file');
    audiowrite(FileName, handles.current_audio, handles.playlist.SampleRate);
    
function export_temp_button_Callback(hObject, eventdata, handles)
        % exporte le graphique temporel courant en JPG
    FileName = uiputfile('*.jpg', 'Save your fantastic temporal chart');
    f=getframe(handles.temp);
    [x,map]=frame2im(f);
    imwrite(x,FileName,'jpg');


function export_freq_button_Callback(hObject, eventdata, handles)
        % exporte le graphique frequentiel courant en JPG
    FileName = uiputfile('*.jpg', 'Save your fantastic temporal chart');
    f=getframe(handles.freq);
    [x,map]=frame2im(f);
    imwrite(x,FileName,'jpg');

% --- Executes on button press in export_poles.
function export_poles_Callback(hObject, eventdata, handles)
        % exporte le graphique de poles zeros courant en JPG
    FileName = uiputfile('*.jpg', 'Save your fantastic PZMAP chart');
    f=getframe(handles.poleszero);
    [x,map]=frame2im(f);
    imwrite(x,FileName,'jpg');

    

% --- Executes on button press in export_bode.
function export_bode_Callback(hObject, eventdata, handles)
        % exporte le graphique de bode courant en JPG
    FileName = uiputfile('*.jpg', 'Save your fantastic Bode chart');
    f=getframe(handles.bode);
    [x,map]=frame2im(f);
    imwrite(x,FileName,'jpg');

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% FONCTIONS INUTILISEES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
function echoslider_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
function freqslider_CreateFcn(hObject, eventdata, handles)
    
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

    
function volumeslider_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function passe_bas_slider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function slider5_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
function timeslider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function original_checkbox_Callback(hObject, eventdata, handles)

function original_freq_checkbox_Callback(hObject, eventdata, handles)


% --- Executes on slider movement.

function highpassslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highpassslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.

% --- Executes during object creation, after setting all properties.
function ordreslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ordreslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



    
